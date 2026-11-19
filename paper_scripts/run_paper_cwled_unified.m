function T = run_paper_cwled_unified()
%RUN_PAPER_CWLED_UNIFIED Four architectures under the digitized Klitzke CW-LED spectrum scaled to 3.09 W/m2 (1000 lx), 0.02 V grid, unified extractor (all finite points,
% J_sc at 0 V, linear V_oc, parabolic MPP); same setup as run_paper_indoor_current_1000lx.m otherwise. Output: results/paper_indoor_cwled_unified/cwled_1000lx.csv
assertCanonicalPath(); warning('off','all');
root = fileparts(fileparts(mfilename('fullpath'))); outDir = fullfile(root,'results','paper_indoor_cwled_unified'); if ~isfolder(outDir), mkdir(outDir); end
arc = load(fullfile(root,'results','paper_stageA','stageA_comparison.mat'),'d1nm','d2nm');
arcFn = @(lam) tmmReflectance(lam,[1.0,1.38,2.35,gapRefractiveIndexReal(lam)],[arc.d1nm,arc.d2nm]);
Bmat = paperRadiativeCoefficients(300);
Tspd = readtable(resolveDataFile('Klitzke2025_CWLED_SPD_digitized.csv'));
spec = [Tspd.lambda_nm, Tspd.relative_irradiance]; spec(:,2) = spec(:,2)*(3.09/spectrumPower(spec)); Pin = spectrumPower(spec);
Vgrid = (0:0.02:1.56)'; rows = {};
for code = ["HJ-0","F-HJ","RBSF-GaP","RHJ"]
  dev = buildPaperArchitecture(char(code), ...
    'BradAbs',Bmat.GaAsPN_cm3ps,'BradGaP',Bmat.GaP_cm3ps,'BradGaNP',Bmat.GaPN_Ege18_cm3ps, ...
    'CnAugerAbs',3e-30,'CpAugerAbs',3e-30,'CnAugerGaP',3e-30,'CpAugerGaP',3e-30,'CnAugerGaNP',3e-30,'CpAugerGaNP',3e-30);
  G = deviceGenAndTransmitML(dev,spec,arcFn);
  gate = indoorSweepLocal(dev,G,Vgrid,Pin);
  f = gate.isFinite(:); pts = [gate.Vs(f), gate.Jterm(f)]; m = pvMetrics(pts,Pin);
  writematrix(pts,fullfile(outDir,sprintf('%s_cwled_1000lx_jv.csv',char(code))));
  rows(end+1,:) = {char(code), Pin, m.Jsc*1e6, m.Voc*1e3, m.FF*100, m.PCE, char(gate.curveClass)}; %#ok<AGROW>
  T = cell2table(rows,'VariableNames',{'Architecture','Pin_Wm2','Jsc_uAcm2','Voc_mV','FF_pct','PCE_pct','Class'});
  writetable(T,fullfile(outDir,'cwled_1000lx.csv')); disp(T(end,:));
end
end

function gate = indoorSweepLocal(dev,G,Vgrid,Pin)
N = dev.N; idxStart = find(abs(Vgrid-0.02)<1e-12,1);
r002 = gummelSolveQFPAnalytic(dev,0.02,G,struct('psiInit',dev.psiEq,'phinInit',zeros(N,1),'phipInit',zeros(N,1))); prev = r002;
for v = [0.01,0.005,0], r0 = gummelSolveQFPAnalytic(dev,v,G,struct('psiInit',prev.psi,'phinInit',prev.phin,'phipInit',prev.phip)); prev = r0; end
states = repmat(r002,numel(Vgrid),1); states(1) = r0; states(idxStart) = r002; prev = r002;
for k = idxStart+1:numel(Vgrid)
  rk = gummelSolveQFPAnalytic(dev,Vgrid(k),G,struct('psiInit',prev.psi,'phinInit',prev.phin,'phipInit',prev.phip));
  states(k) = rk; if rk.isFinite, prev = rk; end
end
gate = indoorAcceptanceGate(dev,G,Vgrid,states,Pin);
end
