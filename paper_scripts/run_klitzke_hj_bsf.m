function T = run_klitzke_hj_bsf(cond, bsf)
%RUN_KLITZKE_HJ_BSF Klitzke GaInP homojunction WITH the back-surface-field layer shown in Klitzke et al. Fig. 1(a)
% (the earlier HJ model, run_klitzke_hj_rhj_production.m, had the p base directly on the rear contact).
% Everything else as in run_klitzke_hj_rhj_production.m (indoor: FD-Jacobian QFP solver, 0.02 V grid to 1.50 V,
% indoorAcceptanceGate, pvMetrics on all finite points) and run_klitzke_table_v112.m 'am15d' (AM1.5D, grid to 1.60 V).
% The BSF composition and doping are not given by Klitzke et al.:
%   bsf = "AlInP"  100 nm p-AlInP, N_A = 2e17 (as the RHJ rear layer of the model)   [default]
%   bsf = "GaInP"  100 nm p+-GaInP, N_A = 1e18 (homo-BSF; weaker barrier, sensitivity)
% cond: "100", "1000" (lx, CW-LED 0.31/3.09 W/m2) or "AM15D".
% Output: results/paper_klitzke_hj_bsf/HJBSF_<bsf>_<cond>.csv (metrics) and _curve.csv
if nargin < 2 || isempty(bsf), bsf = "AlInP"; end
cond = string(cond); bsf = string(bsf);
assertCanonicalPath(); warning('off','all');
root = fileparts(fileparts(mfilename('fullpath'))); out = fullfile(root,'results','paper_klitzke_hj_bsf'); if ~isfolder(out), mkdir(out); end
alOpt = loadAlInPOptical(); arcFn = @(lam) tmmReflectance(lam,[1,1.38,2.35,complexIndex(alOpt,lam)],[85,45]);
exp_ = struct('c100',[11.5 1040 64.2 24.7],'c1000',[114 1240 75.5 34.5],'cAM15D',[14.1e3 1440 86.7 17.6],'cAM15D1000',[14.1e3 1440 86.7 17.6]);   % Klitzke Table II (Jsc in uA/cm2)
L = homojunctionLayers(alOpt);
if bsf == "AlInP", L(4) = layer(100e-7,30,2.35,3.78,11.8,0,2e17,100,10,1e-9,2.5e18,7e18,0,.05);
else, L(4) = layer(100e-7,30,1.90,4.10,11.8,0,1e18,500,30,1e-9,6.5e17,1.5e19,3e4,.05); end
dev = setupDeviceML(L,3,struct('Sn_int',{0,0,0},'Sp_int',{0,0,0}),struct('Sback',1e5));
if cond == "AM15D"
  spec = loadRealAM15D(); Vg = (0:0.02:1.60)';
elseif cond == "AM15D1000"
  % ASTM G173 direct integrates to 900 W/m2; Klitzke et al. measured under AM1.5d calibrated to 1000 W/m2
  spec = loadRealAM15D(); spec(:,2) = spec(:,2)*1000/spectrumPower(spec); Vg = (0:0.02:1.60)';
else
  Tspd = readtable(resolveDataFile('Klitzke2025_CWLED_SPD_digitized.csv')); spec = [Tspd.lambda_nm, Tspd.relative_irradiance];
  lx = str2double(cond); spec(:,2) = spec(:,2)*((0.31*(lx==100)+3.09*(lx==1000))/spectrumPower(spec)); Vg = (0:0.02:1.50)';
end
Pin = spectrumPower(spec); G = deviceGenAndTransmitML(dev,spec,arcFn); t0 = tic;
[gate, states] = indoorSweep(dev,G,Vg,Pin);
fU = gate.isFinite(:); m = pvMetrics([gate.Vs(fU), gate.Jterm(fU)],Pin); e = exp_.("c"+cond);
T = table("HJ+BSF("+bsf+")", cond, m.Jsc*1e6, m.Voc*1e3, m.FF*100, m.PCE, string(gate.curveClass), sum([states.numericalConverged]), numel(Vg), e(1), e(2), e(3), e(4), toc(t0), ...
  'VariableNames',{'Device','Condition','Jsc_uAcm2','Voc_mV','FF_pct','PCE_pct','Class','nConverged','nV','ExpJsc','ExpVoc','ExpFF','ExpPCE','Seconds'});
writetable(T,fullfile(out,sprintf('HJBSF_%s_%s.csv',bsf,cond)));
writetable(table(gate.Vs(:), gate.Jterm(:), fU, 'VariableNames',{'V_V','J_Acm2','Finite'}), fullfile(out,sprintf('HJBSF_%s_%s_curve.csv',bsf,cond)));
disp(T)
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
function L=layer(t,N,Eg,chi,eps_,Nd,Na,muN,muP,tau,Nc,Nv,alphaMax,edgeWidth)
L=struct('t',t,'N',N,'Eg',Eg,'chi',chi,'eps',eps_,'Nd',Nd,'Na',Na,'muN',muN,'muP',muP,'tau',tau,'Nc',Nc,'Nv',Nv,'NcNv',sqrt(Nc*Nv),'Brad',0,'CnAuger',3e-30,'CpAuger',3e-30,'alphaMax',alphaMax,'edgeWidth',edgeWidth,'useOpticalTable',false,'opticalTable',[],'srhCenters',[]);
end
