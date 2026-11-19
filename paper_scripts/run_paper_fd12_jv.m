function S = run_paper_fd12_jv(code, icond)
%RUN_PAPER_FD12_JV Fermi-Dirac (F_1/2) J-V curves by homotopy in the statistics (round-5 review, 2026-11-08).
% At every voltage the F_1/2 solution (consistent degenerate SG, assembleResidualQFPFD12.m with the expm1 Bernoulli function) is
% started from the converged production Boltzmann state at the same voltage, with psi shifted by the difference of the two
% equilibrium potentials. Production physics (Stage-A ARC, SRH + radiative + Auger, nominal lifetime), 0.02 V grid.
% icond: 1 = AM1.5G, 2 = CIE 6500 K / 1000 lx, 3 = CIE 6500 K / 100 lx.
% Both curves are evaluated with the same unified extractor (pvMetrics on all finite points).
% Output: results/paper_fd12_jv/<code>_<cond>.csv (V, J_B, J_FD, flags) and <code>_<cond>_summary.csv.
assertCanonicalPath(); warning('off','all');
root = fileparts(fileparts(mfilename('fullpath'))); outDir = fullfile(root,'results','paper_fd12_jv'); if ~isfolder(outDir), mkdir(outDir); end
condNames = {'AM1.5G','6500K_1000lx','6500K_100lx'};
switch icond
  case 1, spec = loadRealAM15G(); Vg = (0:0.02:1.70)';
  case 2, spec = buildCIEIndoorSpectrum('6500K',1000); Vg = (0:0.02:1.56)';
  case 3, spec = buildCIEIndoorSpectrum('6500K',100); Vg = (0:0.02:1.50)';
end
Pin = spectrumPower(spec);
Bmat = paperRadiativeCoefficients(300);
common = {'BradAbs',Bmat.GaAsPN_cm3ps,'BradGaP',Bmat.GaP_cm3ps,'BradGaNP',Bmat.GaPN_Ege18_cm3ps,'CnAugerAbs',3e-30,'CpAugerAbs',3e-30,'CnAugerGaP',3e-30,'CpAugerGaP',3e-30,'CnAugerGaNP',3e-30,'CpAugerGaNP',3e-30};
arc = load(fullfile(root,'results','paper_stageA','stageA_comparison.mat'),'d1nm','d2nm');
arcFn = @(lam) tmmReflectance(lam,[1.0,1.38,2.35,gapRefractiveIndexReal(lam)],[arc.d1nm,arc.d2nm]);
dev = buildPaperArchitecture(code, common{:}); N = dev.N; Z = zeros(N,1);
psiF = solveEquilibriumFD12(dev); d = dev; d.psiEq = psiF;
[n,p] = carrierDensitiesQFPFD12(d, psiF, Z, Z); d.ni = sqrt(n.*p);
if isfield(d,'srhN1'), d.srhN1(:,1) = d.ni; d.srhP1(:,1) = d.ni; end
G = deviceGenAndTransmitML(dev,spec,arcFn);
% Boltzmann sweep in the order of the production indoor sweep: 0.02 V cold, down to 0 V, then forward
order = [find(abs(Vg-0.02)<1e-12); 1; (find(abs(Vg-0.02)<1e-12)+1:numel(Vg))'];
nV = numel(Vg); JB = nan(nV,1); JF = nan(nV,1); cB = false(nV,1); cF = false(nV,1); fB = false(nV,1); fF = false(nV,1); dJF = nan(nV,1); itF = nan(nV,1); tF = nan(nV,1);
prev = struct('psi',dev.psiEq,'phin',Z,'phip',Z); r02 = [];
for k = order(:)'
  if k == 1 && ~isempty(r02), prev = r02; end
  rB = gummelSolveQFPAnalytic(dev,Vg(k),G,struct('psiInit',prev.psi,'phinInit',prev.phin,'phipInit',prev.phip));
  if abs(Vg(k)-0.02) < 1e-12, r02 = rB; end
  JB(k) = rB.J; cB(k) = rB.numericalConverged; fB(k) = rB.isFinite;
  if rB.isFinite, prev = rB; end
  t0 = tic;
  rF = gummelSolveQFPFD12(d,Vg(k),G,struct('psiInit',rB.psi + (psiF - dev.psiEq),'phinInit',rB.phin,'phipInit',rB.phip,'maxOuter',300));
  JF(k) = rF.J; cF(k) = rF.numericalConverged; fF(k) = rF.isFinite; dJF(k) = rF.deltaJ; itF(k) = rF.iters; tF(k) = toc(t0);
  fprintf('%s %s V=%.2f JB=%.6e JF=%.6e convB=%d convF=%d dJF=%.2e it=%d (%.0f s)\n', code, condNames{icond}, Vg(k), JB(k), JF(k), cB(k), cF(k), dJF(k), itF(k), tF(k));
  if k == 1, prev = r02; end
end
T = table(Vg, JB, JF, cB, cF, fB, fF, dJF, itF, tF, 'VariableNames', {'V_V','JB_Acm2','JFD_Acm2','convB','convFD','finiteB','finiteFD','deltaJ_FD','itersFD','secondsFD'});
writetable(T, fullfile(outDir, sprintf('%s_%s.csv', code, condNames{icond})));
mB = pvMetrics([Vg(fB), JB(fB)], Pin); mF = pvMetrics([Vg(fF), JF(fF)], Pin);
S = table(string(code), string(condNames{icond}), mB.Jsc, mF.Jsc, mB.Voc, mF.Voc, mB.FF, mF.FF, mB.PCE, mF.PCE, sum(cB), sum(cF), nV, max(dJF(fF)), ...
  'VariableNames', {'Architecture','Condition','JscB','JscFD','VocB','VocFD','FFB','FFFD','PCEB','PCEFD','nConvB','nConvFD','nV','maxDeltaJ_FD'});
writetable(S, fullfile(outDir, sprintf('%s_%s_summary.csv', code, condNames{icond})));
disp(S);
end
