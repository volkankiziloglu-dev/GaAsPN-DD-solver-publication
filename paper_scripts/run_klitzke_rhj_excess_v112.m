function run_klitzke_rhj_excess_v112(meshScale)
%RUN_KLITZKE_RHJ_EXCESS_V112 Excess current of the lower-doped Klitzke RHJ with this project's device definition (dd_solver v1.1.2).
% dJ(V) = J_rec,meas - J_rec,model, J_rec = Jsc - J(V), at 100, 200, 300, 500 and 1000 lx (CW-LED scaled linearly from 3.09 W/m2 at
% 1000 lx). Measured curves: data/experimental/klitzke2025/klitzke2025_fig4_iv_curves.csv (LowerRHJ; digitized Fig. 4 of Klitzke et al.,
% copied from Klitzke_Validasyon). Model: RHJ N_D = 2e16, rear N_A = 2e17 (as run_klitzke_table_v112.m), AlInP optical table, 85/45 nm
% ARC, analytic-Jacobian QFP production solver, intrinsic. Measured Jsc = mean over 0.06-0.30 V; model Jsc = J(0).
% Replaces the table computed in Klitzke_Validasyon (klitzke_rhj_light_dependence.m: N_D = 5e17, other optics).
% Output: results/paper_klitzke_rhj_excess_v112/ (model curves, rhj_excess_current.csv, summary).
if nargin < 1, meshScale = 4; end
root = fileparts(fileparts(mfilename('fullpath')));
restoredefaultpath; addpath(fullfile(root,'src'),'-begin'); assertCanonicalPath(); warning('off','all');
out = fullfile(root,'results','paper_klitzke_rhj_excess_v112'); if ~isfolder(out), mkdir(out); end
alOpt = loadAlInPOptical(); arcFn = @(lam) tmmReflectance(lam,[1,1.38,2.35,cidx(alOpt,lam)],[85,45]);
Tspd = readtable(resolveDataFile('Klitzke2025_CWLED_SPD_digitized.csv')); led = [Tspd.lambda_nm, Tspd.relative_irradiance];
Tm = readtable(resolveDataFile(fullfile('experimental','klitzke2025','klitzke2025_fig4_iv_curves.csv')));
L(1)=layer(25e-7,30,2.35,3.78,11.8,8e18,0,100,10,1e-9,2.5e18,7e18,0,.05); L(1).useOpticalTable=true; L(1).opticalTable=alOpt;
L(2)=layer(850e-7,190,1.90,4.10,11.8,2e16,0,500,30,248e-9,6.5e17,1.5e19,3e4,.05); L(2).Brad=1e-10; L(2).srhCenters=struct('tauN',337e-9,'tauP',248e-9,'EtRelEi_eV',0);
L(3)=layer(100e-7,30,2.35,3.78,11.8,0,2e17,100,10,1e-9,2.5e18,7e18,0,.05);
for k = 1:3, L(k).N = max(2,round(L(k).N*meshScale)); end
dev = setupDeviceML(L,2,struct('Sn_int',{0,0},'Sp_int',{0,0}),struct('Sback',1e5)); N = dev.N;
lxs = [100 200 300 500 1000]; Vg = (0:0.01:1.50)'; Vq = (0.70:0.05:1.25)';
dJ = nan(numel(Vq),numel(lxs)); S = nan(numel(lxs),8);
for i = 1:numel(lxs)
  lx = lxs(i); spec = led; spec(:,2) = spec(:,2)*(3.09e-3*lx/spectrumPower(spec)); Pin = spectrumPower(spec);
  G = deviceGenAndTransmitML(dev,spec,arcFn);
  i2 = find(abs(Vg-0.02)<1e-12,1);
  r02 = gummelSolveQFPAnalytic(dev,0.02,G,struct('psiInit',dev.psiEq,'phinInit',zeros(N,1),'phipInit',zeros(N,1))); prev = r02;
  for v = [0.01 0.005 0], r0 = gummelSolveQFPAnalytic(dev,v,G,struct('psiInit',prev.psi,'phinInit',prev.phin,'phipInit',prev.phip)); prev = r0; end
  St = repmat(r02,numel(Vg),1); St(1) = r0; St(i2) = r02; prev = r02;
  for k = i2+1:numel(Vg), rk = gummelSolveQFPAnalytic(dev,Vg(k),G,struct('psiInit',prev.psi,'phinInit',prev.phin,'phipInit',prev.phip)); St(k) = rk; if rk.isFinite, prev = rk; end, end
  g = indoorAcceptanceGate(dev,G,Vg,St,Pin); f = g.isFinite(:); m = pvMetrics([g.Vs(f), g.Jterm(f)],Pin);
  mo = [g.Vs(f), g.Jterm(f)*1e6]; writematrix(mo, fullfile(out,sprintf('rhj_model_%dlx_mesh%d.csv',lx,meshScale)));
  sel = strcmp(Tm.device,'LowerRHJ') & Tm.lux == lx; [v,kk] = sort(Tm.V_V(sel)); j = Tm.J_uAcm2(sel); me = [v, j(kk)];
  [~,iu] = unique(me(:,1)); me = me(iu,:);
  JscM = mean(me(me(:,1)>0.06 & me(:,1)<0.3,2)); JscD = mo(1,2);
  dJ(:,i) = (JscM - interp1(me(:,1),me(:,2),Vq,'linear',NaN)) - (JscD - interp1(mo(:,1),mo(:,2),Vq,'linear',NaN));
  S(i,:) = [lx, Pin, JscM, m.Jsc*1e6, m.Voc*1e3, m.FF*100, m.PCE, sum([St.numericalConverged])];
  fprintf('%4d lx: Jsc meas %.3f model %.3f uA, Voc %.1f mV, FF %.2f, PCE %.3f, %s %d/%d\n',lx,JscM,m.Jsc*1e6,m.Voc*1e3,m.FF*100,m.PCE,char(g.curveClass),S(i,8),numel(Vg));
end
writetable(array2table([Vq dJ],'VariableNames',['V_V', strcat('dJ_',string(lxs))]), fullfile(out,sprintf('rhj_excess_current_mesh%d.csv',meshScale)));
writetable(array2table(S,'VariableNames',{'Lux','Pin_Wm2','JscMeas_uAcm2','Jsc_uAcm2','Voc_mV','FF_pct','PCE_pct','nConverged'}), fullfile(out,sprintf('rhj_model_summary_mesh%d.csv',meshScale)));
disp(array2table([Vq dJ],'VariableNames',['V_V', strcat('dJ_',string(lxs))]));
end
function z=cidx(opt,lam)
E=1239.842/lam; n=interp1(opt.E,opt.n,E,'linear','extrap'); a=interp1(opt.E,opt.alpha,E,'linear',0); k=max(a,0)*(lam*1e-7)/(4*pi); z=n-1i*k;
end
function L=layer(t,N,Eg,chi,eps_,Nd,Na,muN,muP,tau,Nc,Nv,alphaMax,edgeWidth)
L=struct('t',t,'N',N,'Eg',Eg,'chi',chi,'eps',eps_,'Nd',Nd,'Na',Na,'muN',muN,'muP',muP,'tau',tau,'Nc',Nc,'Nv',Nv,'NcNv',sqrt(Nc*Nv),'Brad',0,'CnAuger',3e-30,'CpAuger',3e-30,'alphaMax',alphaMax,'edgeWidth',edgeWidth,'useOpticalTable',false,'opticalTable',[],'srhCenters',[]);
end
