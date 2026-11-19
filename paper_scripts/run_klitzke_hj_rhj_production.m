function T = run_klitzke_hj_rhj_production(meshScale)
%RUN_KLITZKE_HJ_RHJ_PRODUCTION Klitzke GaInP homojunction (HJ) and lower-doped RHJ with ONE solver (round-5 review R1/R5):
% the quasi-Fermi-potential solver with the production residual and a finite-difference Jacobian (gummelSolveQFP; the analytic-Jacobian
% variant used in production does not handle the three SRH centres of the HJ base and fails at the first point), the production indoor sweep (0.02 V grid, 0.02 V cold start,
% continuation to 0 V, forward continuation), indoorAcceptanceGate class and the unified extractor (pvMetrics on all finite points).
% Devices and optics as in run_klitzke2025_indoor_alinp_arc_trial.m (AlInP window optical table, 85/45 nm MgF2/ZnS ARC; HJ with the three
% p-base SRH centres of Klitzke et al. Table I; RHJ N_D = 2e16, rear N_A = 2e17, single centre 337/248 ns, B = 1e-10);
% digitized CW-LED SPD scaled to 0.31 and 3.09 W/m2 (100 and 1000 lx). For reference the HJ is also solved with the earlier carrier-density
% Gummel sweep (sweepJVML, 0.01 V) used for the HJ rows of the validation table.
% Output: results/paper_klitzke_one_solver/klitzke_one_solver_mesh<k>.csv and the J-V curves.
if nargin < 1, meshScale = 1; end
assertCanonicalPath(); warning('off','all');
root = fileparts(fileparts(mfilename('fullpath'))); out = fullfile(root,'results','paper_klitzke_one_solver'); if ~isfolder(out), mkdir(out); end
alOpt = loadAlInPOptical(); arcFn = @(lam) tmmReflectance(lam,[1,1.38,2.35,complexIndex(alOpt,lam)],[85,45]);
Tspd = readtable(resolveDataFile('Klitzke2025_CWLED_SPD_digitized.csv')); specRaw = [Tspd.lambda_nm, Tspd.relative_irradiance];
exp_ = struct('HJ100',[11.5 1040 64.2 24.7],'HJ1000',[114 1240 75.5 34.5],'RHJ100',[11.3 1280 81.1 37.5],'RHJ1000',[112 1350 83.4 41.0]);
rows = {};
for dn = ["HJ","RHJ"]
  if dn == "HJ", L = homojunctionLayers(alOpt); ia = 3; else, L = rhjLayers(2e16,2e17,alOpt); ia = 2; end
  for k = 1:numel(L), L(k).N = max(2,round(L(k).N*meshScale)); end
  dev = setupDeviceML(L,ia,struct('Sn_int',{0,0},'Sp_int',{0,0}),struct('Sback',1e5)); N = dev.N;
  for lx = [100 1000]
    spec = specRaw; spec(:,2) = spec(:,2)*((0.31*(lx==100)+3.09*(lx==1000))/spectrumPower(spec)); Pin = spectrumPower(spec);
    G = deviceGenAndTransmitML(dev,spec,arcFn);
    Vg = (0:0.02:1.50)'; t0 = tic;
    [gate, states] = indoorSweep(dev,G,Vg,Pin);
    fU = gate.isFinite(:); m = pvMetrics([gate.Vs(fU), gate.Jterm(fU)],Pin);
    e = exp_.(sprintf('%s%d',dn,lx));
    rows(end+1,:) = {char(dn), lx, meshScale, 'QFP production', m.Jsc*1e6, m.Voc*1e3, m.FF*100, m.PCE, char(gate.curveClass), sum([states.numericalConverged]), numel(Vg), e(1), e(2), e(3), e(4), toc(t0)}; %#ok<AGROW>
    writetable(table(gate.Vs(:), gate.Jterm(:), fU, 'VariableNames',{'V_V','J_Acm2','Finite'}), fullfile(out,sprintf('%s_%dlx_mesh%g_qfp.csv',dn,lx,meshScale)));
    fprintf('%s %d lx mesh %g QFP: Jsc %.3f uA Voc %.1f FF %.2f PCE %.2f class %s\n', dn, lx, meshScale, m.Jsc*1e6, m.Voc*1e3, m.FF*100, m.PCE, char(gate.curveClass));
    if dn == "HJ" && meshScale == 1
      t0 = tic; Vmax = 1.35*(lx==100) + 1.50*(lx==1000);
      [pts,status] = sweepJVML(dev,G,0,Vmax,.01); mo = pvMetrics(pts,Pin);
      rows(end+1,:) = {char(dn), lx, meshScale, 'Gummel (earlier)', mo.Jsc*1e6, mo.Voc*1e3, mo.FF*100, mo.PCE, sprintf('completed=%d',status.completed), NaN, size(pts,1), e(1), e(2), e(3), e(4), toc(t0)}; %#ok<AGROW>
      fprintf('%s %d lx Gummel (earlier): Jsc %.3f uA Voc %.1f FF %.2f PCE %.2f\n', dn, lx, mo.Jsc*1e6, mo.Voc*1e3, mo.FF*100, mo.PCE);
    end
  end
end
T = cell2table(rows,'VariableNames',{'Device','Lux','MeshScale','Solver','Jsc_uAcm2','Voc_mV','FF_pct','PCE_pct','Class','nConverged','nV','ExpJsc','ExpVoc','ExpFF','ExpPCE','Seconds'});
writetable(T,fullfile(out,sprintf('klitzke_one_solver_mesh%g.csv',meshScale))); disp(T);
end

function [gate, states] = indoorSweep(dev,G,Vgrid,Pin)
N = dev.N; idxStart = find(abs(Vgrid-0.02)<1e-12,1);
r002 = gummelSolveQFP(dev,0.02,G,struct('psiInit',dev.psiEq,'phinInit',zeros(N,1),'phipInit',zeros(N,1)));
prev = r002;
for v = [0.01,0.005,0], r0 = gummelSolveQFP(dev,v,G,struct('psiInit',prev.psi,'phinInit',prev.phin,'phipInit',prev.phip)); prev = r0; end
states = repmat(r002,numel(Vgrid),1); states(1) = r0; states(idxStart) = r002; prev = r002;
for k = idxStart+1:numel(Vgrid)
  rk = gummelSolveQFP(dev,Vgrid(k),G,struct('psiInit',prev.psi,'phinInit',prev.phin,'phipInit',prev.phip));
  states(k) = rk; if rk.isFinite, prev = rk; end
end
gate = indoorAcceptanceGate(dev,G,Vgrid,states,Pin);
end

function z=complexIndex(opt,lam)
E=1239.842/lam; n=interp1(opt.E,opt.n,E,'linear','extrap'); a=interp1(opt.E,opt.alpha,E,'linear',0); k=max(a,0)*(lam*1e-7)/(4*pi); z=n-1i*k;
end
function L=homojunctionLayers(alOpt)
L(1)=layer(25e-7,30,2.35,3.78,11.8,8e18,0,100,10,1e-9,2.5e18,7e18,0,.05); L(1).useOpticalTable=true; L(1).opticalTable=alOpt;
L(2)=layer(150e-7,60,1.90,4.10,11.8,2e18,0,500,30,248e-9,6.5e17,1.5e19,3e4,.05);
L(3)=layer(700e-7,160,1.90,4.10,11.8,0,2e16,500,30,1.66e-9,6.5e17,1.5e19,3e4,.05);
ctr(1)=struct('tauN',1.66e-9,'tauP',6.23e4*1e-9,'EtRelEi_eV',0); ctr(2)=struct('tauN',5.99e-9,'tauP',4.48e3*1e-9,'EtRelEi_eV',0); ctr(3)=struct('tauN',39.5e-9,'tauP',89.6e-9,'EtRelEi_eV',0); L(3).srhCenters=ctr;
end
function L=rhjLayers(ND,rearNA,alOpt)
L(1)=layer(25e-7,30,2.35,3.78,11.8,8e18,0,100,10,1e-9,2.5e18,7e18,0,.05); L(1).useOpticalTable=true; L(1).opticalTable=alOpt;
L(2)=layer(850e-7,190,1.90,4.10,11.8,ND,0,500,30,248e-9,6.5e17,1.5e19,3e4,.05); L(2).Brad=1e-10; L(2).srhCenters=struct('tauN',337e-9,'tauP',248e-9,'EtRelEi_eV',0);
L(3)=layer(100e-7,30,2.35,3.78,11.8,0,rearNA,100,10,1e-9,2.5e18,7e18,0,.05);
end
function L=layer(t,N,Eg,chi,eps_,Nd,Na,muN,muP,tau,Nc,Nv,alphaMax,edgeWidth)
L=struct('t',t,'N',N,'Eg',Eg,'chi',chi,'eps',eps_,'Nd',Nd,'Na',Na,'muN',muN,'muP',muP,'tau',tau,'Nc',Nc,'Nv',Nv,'NcNv',sqrt(Nc*Nv),'Brad',0,'CnAuger',3e-30,'CpAuger',3e-30,'alphaMax',alphaMax,'edgeWidth',edgeWidth,'useOpticalTable',false,'opticalTable',[],'srhCenters',[]);
end
