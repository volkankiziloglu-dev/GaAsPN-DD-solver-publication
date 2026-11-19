%EXAMPLE_AM15G_RBSF_GAP AM1.5G J-V curve of the RBSF-GaP cell (nominal parameters of the paper).
% n+-GaP window / n+-GaAsPN emitter / p-GaAsPN absorber (3% N, Eg = 1.857 eV, tau = 0.6 us) / p+-GaP rear surface field,
% MgF2/ZnS anti-reflection coating (Stage-A thicknesses), SRH + radiative + Auger recombination, 0.02 V grid.
% The metrics are extracted from all finite points of the sweep (unified extractor).
% Reference (paper, Table "outdoor performance"): Jsc = 16.815 mA/cm2, Voc = 1467.1 mV, FF = 90.66 %, PCE = 22.36 %.
% Run from the package root after startup.m. Runtime: under 1 min on a desktop PC (about 25 s in the smoke test).
% Requires the third-party optical tables listed in docs/data_sources.md (data/external/ or GAASPN_EXTERNAL_DATA).
root = fileparts(fileparts(mfilename('fullpath')));
assertCanonicalPath(); warning('off','MATLAB:nearlySingularMatrix'); warning('off','MATLAB:singularMatrix');

arcFile = fullfile(root,'results','paper_stageA','stageA_comparison.mat');
if ~isfile(arcFile), arcFile = fullfile(root,'expected_outputs','results','paper_stageA','stageA_comparison.mat'); end
arc = load(arcFile,'d1nm','d2nm');                                   % ARC thicknesses (nm) of the Stage-A optimisation
arcFn = @(lam) tmmReflectance(lam,[1.0,1.38,2.35,gapRefractiveIndexReal(lam)],[arc.d1nm,arc.d2nm]);
Bmat = paperRadiativeCoefficients(300);                              % van Roosbroeck-Shockley B_rad of GaAsPN, GaP, GaNP

dev = buildPaperArchitecture("RBSF-GaP", BradAbs=Bmat.GaAsPN_cm3ps, BradGaP=Bmat.GaP_cm3ps, BradGaNP=Bmat.GaPN_Ege18_cm3ps, ...
  CnAugerAbs=3e-30, CpAugerAbs=3e-30, CnAugerGaP=3e-30, CpAugerGaP=3e-30, CnAugerGaNP=3e-30, CpAugerGaNP=3e-30);
spec = loadRealAM15G(); Pin = spectrumPower(spec);                    % ASTM G173 AM1.5G, W/m2
G = deviceGenAndTransmitML(dev, spec, arcFn);                        % TMM/Beer-Lambert generation profile
Vgrid = (0:0.02:1.56)';

t0 = tic;
[curve, diagInfo, m] = runAcceptedJVSweep(dev, G, Vgrid, Pin, true);
tRun = toc(t0);

fprintf('\nRBSF-GaP, AM1.5G (%.1f W/m2), %d of %d points accepted by the gate, %.0f s\n', Pin, diagInfo.nAccepted, diagInfo.nAttempted, tRun);
fprintf('  Jsc = %.3f mA/cm2   Voc = %.1f mV   FF = %.2f %%   PCE = %.2f %%\n', 1e3*m.Jsc, 1e3*m.Voc, 100*m.FF, m.PCE);
fprintf('  paper: PCE = 22.36 %%\n');
outDir = fullfile(root,'results','examples'); if ~isfolder(outDir), mkdir(outDir); end
writetable(curve, fullfile(outDir,'RBSF-GaP_AM15G_jv.csv'));
fprintf('  J-V curve written to results/examples/RBSF-GaP_AM15G_jv.csv\n');
