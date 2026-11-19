function out = gummelSolveQFPAnalyticX(dev, V, Gprofile, opts)
% EXPERIMENTAL (experimental_lowlight): copy of gummelSolveQFPAnalytic.m; the only difference is that the residual call goes to
% assembleResidualQFPX.m (flux form selected by dev.fluxModeX). The production code is unchanged.
%GUMMELSOLVEQFPANALYTIC Same Gummel-alternating QFP solver as
% gummelSolveQFP.m, but each per-block Newton substep uses the exact
% analytic tridiagonal Jacobian (assemblePsiJacobianAnalyticQFP.m/
% assemblePhinJacobianAnalyticQFP.m/assemblePhipJacobianAnalyticQFP.m)
% instead of assembleJacobianFD.m's central-difference approximation.
% The analytic Jacobians were cross-validated against the FD Jacobian to
% ~2.5e-10 relative Frobenius norm (both with and without interface
% recombination active) in tests/test_QFP_analytic_jacobian_crosscheck.m
% before this file was written -- see MATHEMATICAL_MODEL.md N.16.21/N.17.
%
% This is a deliberate parallel file, not an in-place edit of
% gummelSolveQFP.m (project convention: never modify a working solver in
% place while introducing an unproven variant). All outer-loop logic
% (Gummel alternation, relaxation, junctionCouplingCorrection, the
% strict physicsPass acceptance criteria, output fields) is IDENTICAL to
% gummelSolveQFP.m -- copied verbatim -- so the two solvers are
% comparable apples-to-apples on speed and on any convergence/acceptance
% differences (which should be none, since the Jacobians are exact
% derivatives of the SAME residual, just computed a different way; only
% floating-point rounding in the FD's central-difference step can cause
% a bit-level difference).
c = constants();
N = dev.N;

if nargin < 4, opts = struct(); end
if ~isfield(opts,'maxOuter') || isempty(opts.maxOuter), opts.maxOuter = 80; end
if ~isfield(opts,'maxInner') || isempty(opts.maxInner), opts.maxInner = 30; end
if ~isfield(opts,'tolOuter') || isempty(opts.tolOuter), opts.tolOuter = 1e-10; end
if ~isfield(opts,'tolInner') || isempty(opts.tolInner), opts.tolInner = 1e-13; end
if ~isfield(opts,'relax') || isempty(opts.relax), opts.relax = 1.0; end
if ~isfield(opts,'verbose') || isempty(opts.verbose), opts.verbose = false; end
if ~isfield(opts,'rowScale') || isempty(opts.rowScale), opts.rowScale = false; end

if isfield(opts,'psiInit') && ~isempty(opts.psiInit)
  psi = opts.psiInit;
else
  psi = dev.psiEq;
end
psi(1) = dev.psiEq(1);
psi(N) = dev.psiEq(N) + V;
if isfield(opts,'phinInit') && ~isempty(opts.phinInit)
  phin = opts.phinInit;
else
  phin = zeros(N,1);
end
if isfield(opts,'phipInit') && ~isempty(opts.phipInit)
  phip = opts.phipInit;
else
  phip = zeros(N,1);
end

if ~isfield(opts,'tolJabs') || isempty(opts.tolJabs), opts.tolJabs = 1e-7; end
if ~isfield(opts,'tolJrel') || isempty(opts.tolJrel), opts.tolJrel = 1e-3; end

junctionList = zeros(0,2);
if isfield(dev,'junctionNodes') && ~isempty(dev.junctionNodes)
  for jk = 1:numel(dev.junctionNodes)
    if dev.Sn_int_j(jk) > 0 && dev.Sp_int_j(jk) > 0
      iE = dev.junctionNodes(jk);
      junctionList(end+1,:) = [iE, iE+1]; %#ok<AGROW>
    end
  end
elseif isfield(dev,'Sn_int') && isfield(dev,'Sp_int') && dev.Sn_int > 0 && dev.Sp_int > 0 && isfield(dev,'N_ETL')
  junctionList = [dev.N_ETL, dev.N_ETL+1];
end

isFiniteState = true;
converged = false;
outerIter = 0;
Gzero = zeros(N,1);
for outerIter = 1:opts.maxOuter
  psiPrev = psi; phinPrev = phin; phipPrev = phip;

  % --- Sub-problem 1: Poisson in psi alone (phin,phip frozen) ---
  [psiSolved, okPsi] = newtonSingleBlockAnalytic(dev, psi, phin, phip, V, Gprofile, 1, opts.maxInner, opts.tolInner, 0.5, opts.rowScale);
  if ~okPsi, isFiniteState = false; break; end
  psi = psiPrev + opts.relax*(psiSolved - psiPrev);
  psi(1) = dev.psiEq(1);
  psi(N) = dev.psiEq(N) + V;

  % --- Sub-problem 2: electron continuity in phin alone (psi,phip frozen) ---
  [phinSolved, okPhin] = newtonSingleBlockAnalytic(dev, psi, phin, phip, V, Gprofile, 2, opts.maxInner, opts.tolInner, 1.0, opts.rowScale);
  if ~okPhin, isFiniteState = false; break; end
  phin = phinPrev + opts.relax*(phinSolved - phinPrev);

  % --- Sub-problem 3: hole continuity in phip alone (psi,phin frozen) ---
  [phipSolved, okPhip] = newtonSingleBlockAnalytic(dev, psi, phin, phip, V, Gprofile, 3, opts.maxInner, opts.tolInner, 1.0, opts.rowScale);
  if ~okPhip, isFiniteState = false; break; end
  phip = phipPrev + opts.relax*(phipSolved - phipPrev);

  % --- Sub-problem 4 (N.16.18): local fully-coupled junction correction,
  % identical to gummelSolveQFP.m (see that file's comment) ---
  if ~isempty(junctionList)
    [phin, phip] = junctionCouplingCorrection(dev, psi, phin, phip, V, Gprofile, junctionList);
  end

  outerStep = max([max(abs(psi-psiPrev)), max(abs(phin-phinPrev)), max(abs(phip-phipPrev))]);
  if opts.verbose
    [~,nV,pV,JnSegV,JpSegV] = assembleResidualQFPX(dev, psi, phin, phip, V, Gprofile);
    JtotV = ifaceCorrectTotalCurrent(dev, JnSegV+JpSegV, nV, pV);
    fprintf('  outer=%d outerStep=%.4e Jmed=%.6e Jspread=%.4e\n', outerIter, outerStep, median(JtotV), max(JtotV)-min(JtotV));
  end
  if outerStep < opts.tolOuter && outerIter > 2
    [~, nC, pC, JnSegC, JpSegC] = assembleResidualQFPX(dev, psi, phin, phip, V, Gprofile);
    JtotC = ifaceCorrectTotalCurrent(dev, JnSegC + JpSegC, nC, pC);   % unchanged unless 'four'/'cross' (dd_solver v1.1.1)
    JtermC = median(JtotC);
    deltaJabsC = max(JtotC) - min(JtotC);
    deltaJrelC = deltaJabsC / max(abs(JtermC), 1e-6);
    if deltaJabsC < opts.tolJabs || deltaJrelC < opts.tolJrel
      converged = true;
      break;
    end
  end
end

[F, n, p, JnSeg, JpSeg] = assembleResidualQFPX(dev, psi, phin, phip, V, Gprofile);
isFiniteState = isFiniteState && all(isfinite(F)) && all(isfinite(psi)) && ...
  all(isfinite(n)) && all(isfinite(p));


JtotSeg = ifaceCorrectTotalCurrent(dev, JnSeg + JpSeg, n, p);   % unchanged unless 'four'/'cross' (dd_solver v1.1.1)
Jterminal = median(JtotSeg);
Jfloor = 1e-12;
deltaJ = (max(JtotSeg) - min(JtotSeg)) / max(max(abs(JtotSeg)), Jfloor);

idxP = 2:N-1;
chargeScale = c.Q * max([max(dev.Nd(idxP)), max(dev.Na(idxP)), max(n(idxP)), max(p(idxP)), 1e10]);
poissonResidual = max(abs(F(3*(idxP-1)+1))) / chargeScale;
poissonConverged = poissonResidual < 1e-6;

currentConservation = checkCurrentConservation(JtotSeg, Jterminal, 1e-7, 1e-3);
continuityResiduals = checkContinuityResiduals(dev, Gprofile, n, p, JnSeg, JpSeg, psi);

chiRel = dev.chi - dev.chiRef;
EFn = c.VT*log(n./dev.Nc) - psi - chiRel;
EFp = -psi - chiRel - dev.Eg - c.VT*log(p./dev.Nv);
qfpConsistency = struct();
qfpConsistency.maxAbsErrN = max(abs(EFn + phin));
qfpConsistency.maxAbsErrP = max(abs(EFp + phip));
qfpConsistency.pass = qfpConsistency.maxAbsErrN < 1e-9 && qfpConsistency.maxAbsErrP < 1e-9;

equilibriumCheck = struct();
equilibriumCheck.EFnSpan = max(EFn)-min(EFn);
equilibriumCheck.EFpSpan = max(EFp)-min(EFp);
equilibriumCheck.massActionMaxRel = max(abs(n.*p-dev.ni.^2)./max(dev.ni.^2,1e-300));
equilibriumLike = abs(V) < 1e-14 && max(abs(Gprofile)) < 1;
if equilibriumLike
  equilibriumCheck.pass = equilibriumCheck.EFnSpan<1e-6 && ...
      equilibriumCheck.EFpSpan<1e-6 && equilibriumCheck.massActionMaxRel<1e-6;
  physicsPass = equilibriumCheck.pass;
else
  equilibriumCheck.pass = false;
  physicsPass = currentConservation.pass && ...
      continuityResiduals.pass;
end

out.psi = psi; out.n = n; out.p = p; out.J = Jterminal;
out.phin = phin; out.phip = phip;
out.JnSeg = JnSeg; out.JpSeg = JpSeg;
out.negFloorHit = 0; out.negFloorTotal = 0; out.negFloorMaxPerIter = 0; out.subFloorPositiveTotal = 0;
out.minNRaw = min(n); out.minPRaw = min(p);
out.currentConservation = currentConservation;
out.continuityResiduals = continuityResiduals;
out.equilibriumCheck = equilibriumCheck;
out.qfpConsistency = qfpConsistency;
out.physicsPass = physicsPass;
% Keep numerical iteration convergence separate from the illumination-
% specific physics gate. Indoor acceptance uses this field together with
% its own low-current residual scaling instead of reusing physicsPass.
out.numericalConverged = converged && isFiniteState;
out.converged = converged && isFiniteState && physicsPass;
out.iters = outerIter;
out.isFinite = isFiniteState;
out.poissonConverged = poissonConverged;
out.poissonResidual = poissonResidual;
out.poissonFailCount = 0;
out.deltaJ = deltaJ;
out.Jseg = JtotSeg;
end

function [xBlock, ok] = newtonSingleBlockAnalytic(dev, psi, phin, phip, V, Gprofile, blockIdx, maxIter, tol, maxStep, rowScale)
% Newton solve for ONE variable block using its exact analytic Jacobian
% (blockIdx in {1=psi,2=phin,3=phip}), all other blocks held fixed at
% the values passed in. Mirrors gummelSolveQFP.m's newtonSingleBlock
% line-search/backtracking structure exactly, differing only in how J
% is built (analytic assembler vs FD).
N = dev.N;
switch blockIdx
  case 1, xBlock = psi;
  case 2, xBlock = phin;
  case 3, xBlock = phip;
end
ok = true;
rowIdx = blockIdx:3:3*N;
for iter = 1:maxIter
  [psiC, phinC, phipC] = injectBlock(psi, phin, phip, blockIdx, xBlock);
  Ffull = assembleResidualQFPX(dev, psiC, phinC, phipC, V, Gprofile);
  Fblock = Ffull(rowIdx);
  if ~all(isfinite(Fblock))
    ok = false; return;
  end
  switch blockIdx
    case 1
      J = assemblePsiJacobianAnalyticQFP(dev, psiC, phinC, phipC);
    case 2
      J = assemblePhinJacobianAnalyticQFP(dev, psiC, phinC, phipC, V, Gzeroed(Gprofile));
    case 3
      J = assemblePhipJacobianAnalyticQFP(dev, psiC, phinC, phipC, V, Gzeroed(Gprofile));
  end
  if rowScale
    % Algebraically neutral row equilibration.  Each residual equation and
    % its Jacobian row receive the same factor; only the floating-point
    % conditioning of the sparse solve and line-search norm changes.
    rowNorm = full(max(abs(J),[],2));
    good = isfinite(rowNorm) & rowNorm>0;
    if any(good), floorNorm=max(min(rowNorm(good))*1e-12,realmin); else, floorNorm=1; end
    scale = 1./max(rowNorm,floorNorm);
    S = spdiags(scale,0,N,N);
    delta = (S*J) \ (-(scale.*Fblock));
    normF0 = norm(scale.*Fblock,2);
  else
    scale = [];
    delta = J \ (-Fblock);
    normF0 = norm(Fblock, 2);
  end
  if ~all(isfinite(delta))
    ok = false; return;
  end
  delta = max(min(delta, maxStep), -maxStep);
  alpha = 1.0;
  bestAlpha = 0; bestNorm = normF0;
  for bt = 1:8
    xTrial = xBlock + alpha*delta;
    [psiT, phinT, phipT] = injectBlock(psi, phin, phip, blockIdx, xTrial);
    Ftrial = assembleResidualQFPX(dev, psiT, phinT, phipT, V, Gprofile);
    Ftb = Ftrial(rowIdx);
    if all(isfinite(Ftb))
      if rowScale, nF=norm(scale.*Ftb,2); else, nF=norm(Ftb,2); end
      if nF < bestNorm
        bestNorm = nF; bestAlpha = alpha;
      end
      if nF < normF0
        break;
      end
    end
    alpha = alpha/2;
  end
  step = bestAlpha*delta;
  xBlock = xBlock + step;
  if max(abs(step)) < tol && iter > 1
    break;
  end
end
end

function [psiO, phinO, phipO] = injectBlock(psi, phin, phip, blockIdx, xBlock)
psiO = psi; phinO = phin; phipO = phip;
switch blockIdx
  case 1, psiO = xBlock;
  case 2, phinO = xBlock;
  case 3, phipO = xBlock;
end
end

function Gz = Gzeroed(Gprofile)
% Same rationale as gummelSolveQFP.m's Gzero: the Jacobian of Gprofile's
% purely-additive contribution to Fphin/Fphip is exactly zero regardless
% of its value, but a large illumination-scale G can swamp double
% precision when probed numerically. The analytic Jacobians never
% differentiate G at all (Gprofile only appears in the residual, not in
% these Jacobian formulas), so this is technically unnecessary here --
% kept only so assemblePhin/PhipJacobianAnalyticQFP's call signature
% matches assembleResidualQFP's for readability; passing the real
% Gprofile would give an IDENTICAL Jacobian.
Gz = zeros(size(Gprofile));
end

function [phin, phip] = junctionCouplingCorrection(dev, psi, phin, phip, V, Gprofile, junctionList)
% Identical to gummelSolveQFP.m's local function of the same name.
for jk = 1:size(junctionList,1)
  iE = junctionList(jk,1); iA = junctionList(jk,2);
  rows = [3*(iE-1)+2, 3*(iA-1)+2, 3*(iE-1)+3, 3*(iA-1)+3];
  u = [phin(iE); phin(iA); phip(iE); phip(iA)];
  resFun = @(uu) localJunctionResidual(dev, psi, phin, phip, V, Gprofile, iE, iA, uu, rows);
  F = resFun(u);
  if ~all(isfinite(F)), continue; end
  for iter = 1:20
    if max(abs(F)) < 1e-16
      break;
    end
    h = 1e-6 * max(abs(u), 1e-6);
    J = zeros(4,4);
    for kk = 1:4
      uP = u; uP(kk) = uP(kk) + h(kk);
      FP = resFun(uP);
      if ~all(isfinite(FP)), J = []; break; end
      J(:,kk) = (FP - F) / h(kk);
    end
    if isempty(J), break; end
    delta = J \ (-F);
    if ~all(isfinite(delta)), break; end
    delta = max(min(delta, 0.5), -0.5);
    uTrial = u + delta;
    FTrial = resFun(uTrial);
    if ~all(isfinite(FTrial))
      break;
    end
    u = uTrial; F = FTrial;
    if max(abs(delta)) < 1e-13
      break;
    end
  end
  phin(iE) = u(1); phin(iA) = u(2); phip(iE) = u(3); phip(iA) = u(4);
end
end

function F = localJunctionResidual(dev, psi, phinBase, phipBase, V, Gprofile, iE, iA, u, rows)
phinT = phinBase; phipT = phipBase;
phinT(iE) = u(1); phinT(iA) = u(2); phipT(iE) = u(3); phipT(iA) = u(4);
Ffull = assembleResidualQFPX(dev, psi, phinT, phipT, V, Gprofile);
F = Ffull(rows);
end
