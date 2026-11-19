%EXAMPLE_INDOOR_1000LX Indoor J-V curve of the RBSF-GaP cell under the CIE 6500 K LED spectrum (LED-B5) at 1000 lx.
% Same device and settings as example_am15g_rbsf_gap.m and the indoor production runs (paper_scripts/run_paper_indoor_production.m):
% 0.02 V grid, reverse continuation to V = 0, indoor low-current acceptance gate (curve class PASS-A/PASS-B/SCREEN/FAIL),
% metrics from all finite points (unified extractor).
% Reference (paper, intrinsic, 6500 K / 1000 lx): PCE = 45.45 %, Jsc = 129.2 uA/cm2, Voc = 1325 mV, class PASS-A.
% Run from the package root after startup.m. Runtime: under 1 min on a desktop PC (about 12 s in the smoke test).
% Requires the third-party optical tables listed in docs/data_sources.md (data/external/ or GAASPN_EXTERNAL_DATA).
root = fileparts(fileparts(mfilename('fullpath')));
assertCanonicalPath(); warning('off','MATLAB:nearlySingularMatrix'); warning('off','MATLAB:singularMatrix');

arcFile = fullfile(root,'results','paper_stageA','stageA_comparison.mat');
if ~isfile(arcFile), arcFile = fullfile(root,'expected_outputs','results','paper_stageA','stageA_comparison.mat'); end
arc = load(arcFile,'d1nm','d2nm');
arcFn = @(lam) tmmReflectance(lam,[1.0,1.38,2.35,gapRefractiveIndexReal(lam)],[arc.d1nm,arc.d2nm]);
Bmat = paperRadiativeCoefficients(300);

dev = buildPaperArchitecture("RBSF-GaP", BradAbs=Bmat.GaAsPN_cm3ps, BradGaP=Bmat.GaP_cm3ps, BradGaNP=Bmat.GaPN_Ege18_cm3ps, ...
  CnAugerAbs=3e-30, CpAugerAbs=3e-30, CnAugerGaP=3e-30, CpAugerGaP=3e-30, CnAugerGaNP=3e-30, CpAugerGaNP=3e-30);
spec = buildCIEIndoorSpectrum('6500K', 1000); Pin = spectrumPower(spec);   % CIE LED-B5 scaled to 1000 lx, W/m2
G = deviceGenAndTransmitML(dev, spec, arcFn);
Vgrid = (0:0.02:1.56)';

t0 = tic;
N = dev.N; idxStart = find(abs(Vgrid-0.02) < 1e-12, 1);
r002 = gummelSolveQFPAnalytic(dev, 0.02, G, struct('psiInit',dev.psiEq,'phinInit',zeros(N,1),'phipInit',zeros(N,1)));
prev = r002;
for v = [0.01, 0.005, 0]                                               % reverse continuation to short circuit
  r0 = gummelSolveQFPAnalytic(dev, v, G, struct('psiInit',prev.psi,'phinInit',prev.phin,'phipInit',prev.phip)); prev = r0;
end
states = repmat(r002, numel(Vgrid), 1); states(1) = r0; states(idxStart) = r002; prev = r002;
for k = idxStart+1:numel(Vgrid)                                        % forward sweep
  rk = gummelSolveQFPAnalytic(dev, Vgrid(k), G, struct('psiInit',prev.psi,'phinInit',prev.phin,'phipInit',prev.phip));
  states(k) = rk; if rk.isFinite, prev = rk; end
end
gate = indoorAcceptanceGate(dev, G, Vgrid, states, Pin);
f = gate.isFinite(:);
m = pvMetrics([gate.Vs(f), gate.Jterm(f)], Pin);
tRun = toc(t0);

fprintf('\nRBSF-GaP, CIE 6500 K, 1000 lx (%.3f W/m2), curve class %s, %d finite points, %.0f s\n', Pin, char(gate.curveClass), nnz(f), tRun);
fprintf('  Jsc = %.2f uA/cm2   Voc = %.1f mV   FF = %.2f %%   PCE = %.2f %%\n', 1e6*m.Jsc, 1e3*m.Voc, 100*m.FF, m.PCE);
fprintf('  paper: PCE = 45.45 %%\n');
outDir = fullfile(root,'results','examples'); if ~isfolder(outDir), mkdir(outDir); end
writetable(table(gate.Vs(f), gate.Jterm(f), double(gate.pointPass(f)), 'VariableNames', {'V_V','J_Acm2','Accepted'}), ...
  fullfile(outDir,'RBSF-GaP_6500K_1000lx_jv.csv'));
fprintf('  J-V curve written to results/examples/RBSF-GaP_6500K_1000lx_jv.csv\n');
