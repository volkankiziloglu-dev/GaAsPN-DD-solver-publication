function T = run_klitzke_table_v112(parts)
%RUN_KLITZKE_TABLE_V112 Missing rows for the Klitzke validation table, rebuilt with this project's device definition
% (dd_solver v1.1.2; devices as in run_klitzke_hj_rhj_production.m: AlInP optical table, 85/45 nm MgF2/ZnS ARC, HJ with the
% three p-base SRH centres, RHJ N_D = 2e16, rear N_A = 2e17, single centre 337/248 ns, B = 1e-10).
%   'am15d'  HJ and RHJ under AM1.5D (1000 W/m2), mesh 1, gummelSolveQFP (FD Jacobian) sweep 0..1.60 V
%   'mesh'   RHJ at 100 and 1000 lx on mesh 1, 2 and 4 with the analytic-Jacobian production solver (single SRH centre)
% The indoor mesh-1 rows (HJ and RHJ, 100/1000 lx) are those of results/paper_klitzke_one_solver/klitzke_one_solver_mesh1.csv;
% the RHJ 1000 lx row was reproduced on 2026-11-13 (PCE 42.884 %).
% Output: results/paper_klitzke_table_v112/*.csv
if nargin < 1, parts = {'am15d','mesh'}; end
root = fileparts(fileparts(mfilename('fullpath')));
restoredefaultpath; addpath(fullfile(root,'src'),'-begin'); assertCanonicalPath(); warning('off','all');
out = fullfile(root,'results','paper_klitzke_table_v112'); if ~isfolder(out), mkdir(out); end
alOpt = loadAlInPOptical(); arcFn = @(lam) tmmReflectance(lam,[1,1.38,2.35,cidx(alOpt,lam)],[85,45]);
Tspd = readtable(resolveDataFile('Klitzke2025_CWLED_SPD_digitized.csv')); led = [Tspd.lambda_nm, Tspd.relative_irradiance];
rows = {};
if any(strcmp(parts,'am15d'))
  spec = loadRealAM15D(); Pin = spectrumPower(spec);
  for dn = ["HJ","RHJ"]
    dev = makeDev(dn,1,alOpt); G = deviceGenAndTransmitML(dev,spec,arcFn); t0 = tic;
    [m,cls,nc,nV,Vs,J,f] = sweepQFP(dev,G,(0:0.02:1.60)',Pin,@gummelSolveQFP);
    rows(end+1,:) = {char(dn),'AM1.5D',1,'QFP FD',m.Jsc*1e3,m.Voc*1e3,m.FF*100,m.PCE,cls,nc,nV,toc(t0)}; %#ok<AGROW>
    writetable(table(Vs(:),J(:),f,'VariableNames',{'V_V','J_Acm2','Finite'}),fullfile(out,sprintf('%s_AM15D_mesh1.csv',dn)));
    fprintf('%s AM1.5D: Jsc %.3f mA Voc %.1f FF %.2f PCE %.3f %s %d/%d\n',dn,m.Jsc*1e3,m.Voc*1e3,m.FF*100,m.PCE,cls,nc,nV);
    save_(rows,out);
  end
end
if any(strcmp(parts,'am15d_mesh'))
  % HJ under AM1.5D on mesh 2 with the FD-Jacobian QFP solver (mesh row of the Supplement); compared with the mesh-1 row of 'am15d'.
  % The analytic-Jacobian solver returns no finite point for this device under AM1.5D (tried 2026-11-13).
  spec = loadRealAM15D(); Pin = spectrumPower(spec);
  for ms = 2
    dev = makeDev("HJ",ms,alOpt); G = deviceGenAndTransmitML(dev,spec,arcFn); t0 = tic;
    [m,cls,nc,nV,Vs,J,f] = sweepQFP(dev,G,(0:0.02:1.60)',Pin,@gummelSolveQFP);
    rows(end+1,:) = {'HJ','AM1.5D',ms,'QFP FD',m.Jsc*1e3,m.Voc*1e3,m.FF*100,m.PCE,cls,nc,nV,toc(t0)}; %#ok<AGROW>
    writetable(table(Vs(:),J(:),f,'VariableNames',{'V_V','J_Acm2','Finite'}),fullfile(out,sprintf('HJ_AM15D_mesh%d.csv',ms)));
    fprintf('HJ AM1.5D mesh %d: Jsc %.4f mA Voc %.2f FF %.3f PCE %.4f %s %d/%d\n',ms,m.Jsc*1e3,m.Voc*1e3,m.FF*100,m.PCE,cls,nc,nV);
  end
  T = cell2table(rows,'VariableNames',{'Device','Condition','MeshScale','Solver','Jsc_mAcm2','Voc_mV','FF_pct','PCE_pct','Class','nConverged','nV','Seconds'});
  writetable(T,fullfile(out,'klitzke_table_v112_am15d_mesh.csv')); disp(T); return
end
if any(strcmp(parts,'mesh'))
  for ms = [1 2 4]
    dev = makeDev("RHJ",ms,alOpt);
    for lx = [1000 100]
      spec = led; spec(:,2) = spec(:,2)*((0.31*(lx==100)+3.09*(lx==1000))/spectrumPower(spec)); Pin = spectrumPower(spec);
      G = deviceGenAndTransmitML(dev,spec,arcFn); t0 = tic;
      [m,cls,nc,nV,Vs,J,f] = sweepQFP(dev,G,(0:0.02:1.56)',Pin,@gummelSolveQFPAnalytic);
      rows(end+1,:) = {'RHJ',sprintf('%d lx',lx),ms,'QFP analytic',m.Jsc*1e3,m.Voc*1e3,m.FF*100,m.PCE,cls,nc,nV,toc(t0)}; %#ok<AGROW>
      writetable(table(Vs(:),J(:),f,'VariableNames',{'V_V','J_Acm2','Finite'}),fullfile(out,sprintf('RHJ_%dlx_mesh%d.csv',lx,ms)));
      fprintf('RHJ %d lx mesh %d: Jsc %.4f mA Voc %.1f FF %.2f PCE %.3f %s %d/%d\n',lx,ms,m.Jsc*1e3,m.Voc*1e3,m.FF*100,m.PCE,cls,nc,nV);
      save_(rows,out);
    end
  end
end
T = save_(rows,out); disp(T)
end

function T = save_(rows,out)
T = cell2table(rows,'VariableNames',{'Device','Condition','MeshScale','Solver','Jsc_mAcm2','Voc_mV','FF_pct','PCE_pct','Class','nConverged','nV','Seconds'});
writetable(T,fullfile(out,'klitzke_table_v112_extra.csv'));
end

function [m,cls,nc,nV,Vs,J,f] = sweepQFP(dev,G,Vg,Pin,solver)
N = dev.N; i2 = find(abs(Vg-0.02)<1e-12,1);
r02 = solver(dev,0.02,G,struct('psiInit',dev.psiEq,'phinInit',zeros(N,1),'phipInit',zeros(N,1))); prev = r02;
for v = [0.01 0.005 0], r0 = solver(dev,v,G,struct('psiInit',prev.psi,'phinInit',prev.phin,'phipInit',prev.phip)); prev = r0; end
S = repmat(r02,numel(Vg),1); S(1) = r0; S(i2) = r02; prev = r02;
for k = i2+1:numel(Vg), rk = solver(dev,Vg(k),G,struct('psiInit',prev.psi,'phinInit',prev.phin,'phipInit',prev.phip)); S(k) = rk; if rk.isFinite, prev = rk; end, end
g = indoorAcceptanceGate(dev,G,Vg,S,Pin); f = g.isFinite(:); Vs = g.Vs; J = g.Jterm;
m = pvMetrics([Vs(f), J(f)],Pin); cls = char(g.curveClass); nc = sum([S.numericalConverged]); nV = numel(Vg);
end

function dev = makeDev(dn,ms,alOpt)
if dn == "HJ"
  L(1)=layer(25e-7,30,2.35,3.78,11.8,8e18,0,100,10,1e-9,2.5e18,7e18,0,.05); L(1).useOpticalTable=true; L(1).opticalTable=alOpt;
  L(2)=layer(150e-7,60,1.90,4.10,11.8,2e18,0,500,30,248e-9,6.5e17,1.5e19,3e4,.05);
  L(3)=layer(700e-7,160,1.90,4.10,11.8,0,2e16,500,30,1.66e-9,6.5e17,1.5e19,3e4,.05);
  L(3).srhCenters=struct('tauN',{1.66e-9,5.99e-9,39.5e-9},'tauP',{6.23e-5,4.48e-6,89.6e-9},'EtRelEi_eV',{0,0,0}); ia = 3;
else
  L(1)=layer(25e-7,30,2.35,3.78,11.8,8e18,0,100,10,1e-9,2.5e18,7e18,0,.05); L(1).useOpticalTable=true; L(1).opticalTable=alOpt;
  L(2)=layer(850e-7,190,1.90,4.10,11.8,2e16,0,500,30,248e-9,6.5e17,1.5e19,3e4,.05); L(2).Brad=1e-10; L(2).srhCenters=struct('tauN',337e-9,'tauP',248e-9,'EtRelEi_eV',0);
  L(3)=layer(100e-7,30,2.35,3.78,11.8,0,2e17,100,10,1e-9,2.5e18,7e18,0,.05); ia = 2;
end
for k = 1:numel(L), L(k).N = max(2,round(L(k).N*ms)); end
dev = setupDeviceML(L,ia,struct('Sn_int',{0,0},'Sp_int',{0,0}),struct('Sback',1e5));
end

function z=cidx(opt,lam)
E=1239.842/lam; n=interp1(opt.E,opt.n,E,'linear','extrap'); a=interp1(opt.E,opt.alpha,E,'linear',0); k=max(a,0)*(lam*1e-7)/(4*pi); z=n-1i*k;
end
function L=layer(t,N,Eg,chi,eps_,Nd,Na,muN,muP,tau,Nc,Nv,alphaMax,edgeWidth)
L=struct('t',t,'N',N,'Eg',Eg,'chi',chi,'eps',eps_,'Nd',Nd,'Na',Na,'muN',muN,'muP',muP,'tau',tau,'Nc',Nc,'Nv',Nv,'NcNv',sqrt(Nc*Nv),'Brad',0,'CnAuger',3e-30,'CpAuger',3e-30,'alphaMax',alphaMax,'edgeWidth',edgeWidth,'useOpticalTable',false,'opticalTable',[],'srhCenters',[]);
end
