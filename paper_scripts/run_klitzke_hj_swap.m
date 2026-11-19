function T = run_klitzke_hj_swap(variants, luxes)
%RUN_KLITZKE_HJ_SWAP Which part of the RHJ advantage over the homojunction of Klitzke et al. comes from the lifetime and
% which from the missing rear barrier? Devices as in run_klitzke_table_v112.m (mesh 1, CW-LED 0.31/3.09 W/m2,
% 85/45 nm MgF2/ZnS ARC, S_back = 1e5 cm/s), production analytic-Jacobian QFP solver and indoor acceptance gate.
%   'HJ'        published homojunction (p base with the three SRH centres of the DH samples)
%   'HJ_tau'    p base gets the recombination of the RHJ absorber: one midgap centre with minority (electron) lifetime
%               248 ns and majority (hole) lifetime 337 ns, B = 1e-10 cm3/s (roles mapped from the n-type absorber)
%   'HJ_BSF'    published base + 100 nm p-AlInP rear layer (N_A = 2e17, as the RHJ rear layer)
%   'HJ_tau_BSF' both changes
%   'RHJ'       lower-doped RHJ (reference)
% Solver: analytic-Jacobian QFP; FD-Jacobian QFP (as run_klitzke_hj_rhj_production.m) when the analytic one gives no finite point.
% Output: results/paper_klitzke_hj_swap/klitzke_hj_swap_<variants>.csv
if nargin < 1 || isempty(variants), variants = ["HJ","HJ_tau","HJ_BSF","HJ_tau_BSF","RHJ"]; end
if nargin < 2 || isempty(luxes), luxes = [1000 100]; end
root = fileparts(fileparts(mfilename('fullpath')));
assertCanonicalPath(); warning('off','all');
out = fullfile(root,'results','paper_klitzke_hj_swap'); if ~isfolder(out), mkdir(out); end
alOpt = loadAlInPOptical(); arcFn = @(lam) tmmReflectance(lam,[1,1.38,2.35,cidx(alOpt,lam)],[85,45]);
Tspd = readtable(resolveDataFile('Klitzke2025_CWLED_SPD_digitized.csv')); led = [Tspd.lambda_nm, Tspd.relative_irradiance];
rows = {};
for vn = string(variants)
  dev = makeDev(vn,alOpt);
  for lx = luxes
    spec = led; spec(:,2) = spec(:,2)*((0.31*(lx==100)+3.09*(lx==1000))/spectrumPower(spec)); Pin = spectrumPower(spec);
    G = deviceGenAndTransmitML(dev,spec,arcFn); t0 = tic;
    slv = 'analytic';
    try, [m,cls,nc,nV] = sweepQFP(dev,G,(0:0.02:1.56)',Pin,@gummelSolveQFPAnalytic); catch
      slv = 'FD'; [m,cls,nc,nV] = sweepQFP(dev,G,(0:0.02:1.56)',Pin,@gummelSolveQFP); end
    rows(end+1,:) = {char(vn),slv,lx,m.Jsc*1e6,m.Voc*1e3,m.FF*100,m.PCE,cls,nc,nV,toc(t0)}; %#ok<AGROW>
    fprintf('%-11s %4d lx: Jsc %.3f uA Voc %.1f mV FF %.2f PCE %.3f %s %d/%d (%.0f s)\n',vn,lx,m.Jsc*1e6,m.Voc*1e3,m.FF*100,m.PCE,cls,nc,nV,toc(t0));
    T = cell2table(rows,'VariableNames',{'Device','Solver','Lux','Jsc_uAcm2','Voc_mV','FF_pct','PCE_pct','Class','nConverged','nV','Seconds'});
    writetable(T,fullfile(out,sprintf('klitzke_hj_swap_%s.csv',strjoin(string(variants),'_'))));
  end
end
disp(T)
end

function [m,cls,nc,nV] = sweepQFP(dev,G,Vg,Pin,solver)
N = dev.N; i2 = find(abs(Vg-0.02)<1e-12,1);
r02 = solver(dev,0.02,G,struct('psiInit',dev.psiEq,'phinInit',zeros(N,1),'phipInit',zeros(N,1))); prev = r02;
for v = [0.01 0.005 0], r0 = solver(dev,v,G,struct('psiInit',prev.psi,'phinInit',prev.phin,'phipInit',prev.phip)); prev = r0; end
S = repmat(r02,numel(Vg),1); S(1) = r0; S(i2) = r02; prev = r02;
for k = i2+1:numel(Vg), rk = solver(dev,Vg(k),G,struct('psiInit',prev.psi,'phinInit',prev.phin,'phipInit',prev.phip)); S(k) = rk; if rk.isFinite, prev = rk; end, end
g = indoorAcceptanceGate(dev,G,Vg,S,Pin); f = g.isFinite(:);
m = pvMetrics([g.Vs(f), g.Jterm(f)],Pin); cls = char(g.curveClass); nc = sum([S.numericalConverged]); nV = numel(Vg);
end

function dev = makeDev(vn,alOpt)
win = layer(25e-7,30,2.35,3.78,11.8,8e18,0,100,10,1e-9,2.5e18,7e18,0,.05); win.useOpticalTable=true; win.opticalTable=alOpt;
rear = layer(100e-7,30,2.35,3.78,11.8,0,2e17,100,10,1e-9,2.5e18,7e18,0,.05);
if vn == "RHJ"
  L(1) = win;
  L(2) = layer(850e-7,190,1.90,4.10,11.8,2e16,0,500,30,248e-9,6.5e17,1.5e19,3e4,.05); L(2).Brad=1e-10; L(2).srhCenters=struct('tauN',337e-9,'tauP',248e-9,'EtRelEi_eV',0);
  L(3) = rear; ia = 2;
else
  L(1) = win;
  L(2) = layer(150e-7,60,1.90,4.10,11.8,2e18,0,500,30,248e-9,6.5e17,1.5e19,3e4,.05);
  L(3) = layer(700e-7,160,1.90,4.10,11.8,0,2e16,500,30,1.66e-9,6.5e17,1.5e19,3e4,.05);
  if contains(vn,"tau")
    L(3).tau = 248e-9; L(3).Brad = 1e-10; L(3).srhCenters = struct('tauN',248e-9,'tauP',337e-9,'EtRelEi_eV',0);
  else
    L(3).srhCenters = struct('tauN',{1.66e-9,5.99e-9,39.5e-9},'tauP',{6.23e-5,4.48e-6,89.6e-9},'EtRelEi_eV',{0,0,0});
  end
  if contains(vn,"BSF"), L(4) = rear; end
  ia = 3;
end
dev = setupDeviceML(L,ia,struct('Sn_int',repmat({0},1,numel(L)-1),'Sp_int',repmat({0},1,numel(L)-1)),struct('Sback',1e5));
end

function z=cidx(opt,lam)
E=1239.842/lam; n=interp1(opt.E,opt.n,E,'linear','extrap'); a=interp1(opt.E,opt.alpha,E,'linear',0); k=max(a,0)*(lam*1e-7)/(4*pi); z=n-1i*k;
end
function L=layer(t,N,Eg,chi,eps_,Nd,Na,muN,muP,tau,Nc,Nv,alphaMax,edgeWidth)
L=struct('t',t,'N',N,'Eg',Eg,'chi',chi,'eps',eps_,'Nd',Nd,'Na',Na,'muN',muN,'muP',muP,'tau',tau,'Nc',Nc,'Nv',Nv,'NcNv',sqrt(Nc*Nv),'Brad',0,'CnAuger',3e-30,'CpAuger',3e-30,'alphaMax',alphaMax,'edgeWidth',edgeWidth,'useOpticalTable',false,'opticalTable',[],'srhCenters',[]);
end
