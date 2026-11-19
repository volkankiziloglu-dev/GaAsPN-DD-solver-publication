function summary = run_paper_indoor_current_1000lx()
%RUN_PAPER_INDOOR_CURRENT_1000LX Current-code indoor publication rerun.
% Four main GaAsPN architectures are solved at 1000 lx for the four CIE
% spectra and under the digitized Klitzke CW-LED spectrum scaled to the
% reported 1000-lx incident power (3.09 W/m2). Existing production outputs
% are not overwritten. Point, integral-CV, and native gates are evaluated
% from the same solved states.

assertCanonicalPath();
warning('off','MATLAB:nearlySingularMatrix');
warning('off','MATLAB:singularMatrix');

root = fileparts(fileparts(mfilename('fullpath')));
out = fullfile(root,'results','paper_indoor_current_1000lx');
if ~isfolder(out), mkdir(out); end

arc = load(fullfile(root,'results','paper_stageA','stageA_comparison.mat'),'d1nm','d2nm');
arcFn = @(lam)tmmReflectance(lam,[1,1.38,2.35,gapRefractiveIndexReal(lam)], ...
  [arc.d1nm,arc.d2nm]);
Bmat = paperRadiativeCoefficients(300);

codes = ["HJ-0","F-HJ","RBSF-GaP","RHJ"];
conditions = ["2700K","4000K","6500K","4000F","Klitzke-CWLED"];
dV = 0.01;
V = (0:dV:1.56)';

rows = cell(numel(codes)*numel(conditions),31);
rr = 0;
for ic = 1:numel(conditions)
  [spec,Pin] = makeSpectrum(root,conditions(ic));
  for ia = 1:numel(codes)
    rr = rr + 1;
    fprintf('\n[%d/%d] %s | %s | Pin %.6g W/m2\n', ...
      rr,numel(codes)*numel(conditions),codes(ia),conditions(ic),Pin);

    [dev,meta] = buildPaperArchitecture(codes(ia), ...
      BradAbs=Bmat.GaAsPN_cm3ps,BradGaP=Bmat.GaP_cm3ps, ...
      BradGaNP=Bmat.GaPN_Ege18_cm3ps, ...
      CnAugerAbs=3e-30,CpAugerAbs=3e-30, ...
      CnAugerGaP=3e-30,CpAugerGaP=3e-30, ...
      CnAugerGaNP=3e-30,CpAugerGaNP=3e-30, ...
      HomoRearContactNodes=20);
    G = deviceGenAndTransmitML(dev,spec,arcFn);
    states = sweepStates(dev,G,V);

    gatePoint = indoorAcceptanceGate(dev,G,V,states,Pin,"point");
    gateCV = indoorAcceptanceGate(dev,G,V,states,Pin,"cv");
    gateNative = indoorAcceptanceGate(dev,G,V,states,Pin,"native");
    crit = gateCV.region ~= "deepforward";

    rrPoint = accepted(gatePoint,crit);
    rrCV = accepted(gateCV,crit);
    rrNative = accepted(gateNative,crit);
    rows(rr,:) = {codes(ia),conditions(ic),1000,Pin,dV,dev.N,20, ...
      gateCV.Jsc*1e6,gateCV.Voc*1e3,gateCV.FF*100,gateCV.PCE, ...
      gatePoint.curveClass,rrPoint,sum(crit),gatePoint.scPass,gatePoint.mppPass,gatePoint.vocPass, ...
      gateCV.curveClass,rrCV,gateCV.scPass,gateCV.mppPass,gateCV.vocPass, ...
      gateNative.curveClass,rrNative,gateNative.scPass,gateNative.mppPass,gateNative.vocPass, ...
      safeMax(gateCV.rJ(crit)),safeMax(gateCV.rFB(crit)), ...
      safeMax(gateCV.rnGlobal(crit)),safeMax(gateCV.rpGlobal(crit))};

    C = table(V,gateCV.Jterm,gateCV.region,gatePoint.pointPass,gateCV.pointPass, ...
      gateNative.pointPass,gateCV.rJ,gateCV.rFB,gateCV.rnGlobal,gateCV.rpGlobal, ...
      'VariableNames',{'V_V','J_Acm2','Region','PointPass_point','PointPass_cv', ...
      'PointPass_native','rJ','rFB','rnCV','rpCV'});
    writetable(C,fullfile(out,sprintf('%s_%s_1000lx.csv',codes(ia),conditions(ic))));
    T = makeTable(rows(1:rr,:));
    writetable(T,fullfile(out,'current_1000lx_checkpoint.csv'));
    save(fullfile(out,'current_1000lx_checkpoint.mat'),'T','dV','meta','-v7.3');
    fprintf('  Jsc %.4f uA/cm2, Voc %.2f mV, FF %.3f%%, PCE %.3f%%\n', ...
      gateCV.Jsc*1e6,gateCV.Voc*1e3,gateCV.FF*100,gateCV.PCE);
    fprintf('  point %s %d/%d | CV %s %d/%d | native %s %d/%d\n', ...
      gatePoint.curveClass,rrPoint,sum(crit),gateCV.curveClass,rrCV,sum(crit), ...
      gateNative.curveClass,rrNative,sum(crit));
  end
end

summary = makeTable(rows);
writetable(summary,fullfile(out,'current_1000lx_summary.csv'));
save(fullfile(out,'current_1000lx_summary.mat'),'summary','dV','-v7.3');
disp(summary);
end

function [spec,Pin] = makeSpectrum(root,condition)
if condition == "Klitzke-CWLED"
  T = readtable(resolveDataFile('Klitzke2025_CWLED_SPD_digitized.csv'));
  spec = [T.lambda_nm,T.relative_irradiance];
  Pin = 3.09;
  spec(:,2) = spec(:,2)*(Pin/spectrumPower(spec));
else
  spec = buildCIEIndoorSpectrum(condition,1000);
  Pin = spectrumPower(spec);
end
end

function states = sweepStates(dev,G,V)
N = dev.N;
k02 = find(abs(V-0.02)<1e-12,1);
seed = struct('psi',dev.psiEq,'phin',zeros(N,1),'phip',zeros(N,1));
r02 = gummelSolveQFPAnalytic(dev,0.02,G, ...
  struct('psiInit',seed.psi,'phinInit',seed.phin,'phipInit',seed.phip));
prev = r02;
for v = [0.01,0.005,0]
  r0 = gummelSolveQFPAnalytic(dev,v,G, ...
    struct('psiInit',prev.psi,'phinInit',prev.phin,'phipInit',prev.phip));
  prev = r0;
end
states = repmat(r02,numel(V),1);
states(1) = r0;
states(k02) = r02;
prev = r02;
for k = k02+1:numel(V)
  r = gummelSolveQFPAnalytic(dev,V(k),G, ...
    struct('psiInit',prev.psi,'phinInit',prev.phin,'phipInit',prev.phip));
  states(k) = r;
  if r.isFinite, prev = r; end
end
end

function n = accepted(gate,crit)
n = sum(gate.pointPass & crit);
end

function v = safeMax(x)
x = x(isfinite(x));
if isempty(x), v = NaN; else, v = max(x); end
end

function T = makeTable(rows)
T = cell2table(rows,'VariableNames',{ ...
  'Architecture','Spectrum','Lux','Pin_Wm2','VoltageStep_V','DeviceNodes','HomoRearContactNodes', ...
  'Jsc_uAcm2','Voc_mV','FF_pct','PCE_pct', ...
  'PointClass','PointAccepted','CriticalPoints','PointSC','PointMPP','PointVoc', ...
  'CVClass','CVAccepted','CVSC','CVMPP','CVVoc', ...
  'NativeClass','NativeAccepted','NativeSC','NativeMPP','NativeVoc', ...
  'MaxrJ','MaxrFB','MaxrnCV','MaxrpCV'});
end
