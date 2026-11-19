function out = gummelSolveQFPFD12(dev, V, Gprofile, opts)
%GUMMELSOLVEQFPFD12 (2026-11-07) F_1/2 statistics with the consistent degenerate SG flux (assembleResidualQFPFD12.m); otherwise a copy of gummelSolveQFPFermiDiracDegenSG.m calling the consistent degenerate SG residual assembleResidualQFPFD12.m. Original header: Same as gummelSolveQFPFermiDirac.m, but
% using assembleResidualQFPFD12.m (degeneracy-factor-
% corrected SG diffusion coefficient) instead of the plain-mu*VT version.
% See that file's header and MATHEMATICAL_MODEL.md N.20 for why this is
% an APPROXIMATE (not exact) fix.
c = constants();
N = dev.N;
gamma = 0.27;

if nargin < 4, opts = struct(); end
if ~isfield(opts,'maxOuter') || isempty(opts.maxOuter), opts.maxOuter = 600; end
if ~isfield(opts,'maxInner') || isempty(opts.maxInner), opts.maxInner = 30; end
if ~isfield(opts,'tolOuter') || isempty(opts.tolOuter), opts.tolOuter = 1e-10; end
if ~isfield(opts,'tolInner') || isempty(opts.tolInner), opts.tolInner = 1e-13; end
if ~isfield(opts,'relax') || isempty(opts.relax), opts.relax = 1.0; end
if ~isfield(opts,'verbose') || isempty(opts.verbose), opts.verbose = false; end

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

  resFunPsi = @(psiVec) assembleResidualQFPFD12(dev, psiVec, phin, phip, V, Gprofile);
  jacFunPsi = @(psiVec) assembleResidualQFPFD12(dev, psiVec, phin, phip, V, Gzero);
  [psiSolved, okPsi] = newtonSingleBlock(resFunPsi, jacFunPsi, psi, N, 1, opts.maxInner, opts.tolInner, 0.5);
  if ~okPsi, isFiniteState = false; break; end
  psi = psiPrev + opts.relax*(psiSolved - psiPrev);
  psi(1) = dev.psiEq(1);
  psi(N) = dev.psiEq(N) + V;

  resFunPhin = @(phinVec) assembleResidualQFPFD12(dev, psi, phinVec, phip, V, Gprofile);
  jacFunPhin = @(phinVec) assembleResidualQFPFD12(dev, psi, phinVec, phip, V, Gzero);
  [phinSolved, okPhin] = newtonSingleBlock(resFunPhin, jacFunPhin, phin, N, 2, opts.maxInner, opts.tolInner, 1.0);
  if ~okPhin, isFiniteState = false; break; end
  phin = phinPrev + opts.relax*(phinSolved - phinPrev);

  resFunPhip = @(phipVec) assembleResidualQFPFD12(dev, psi, phin, phipVec, V, Gprofile);
  jacFunPhip = @(phipVec) assembleResidualQFPFD12(dev, psi, phin, phipVec, V, Gzero);
  [phipSolved, okPhip] = newtonSingleBlock(resFunPhip, jacFunPhip, phip, N, 3, opts.maxInner, opts.tolInner, 1.0);
  if ~okPhip, isFiniteState = false; break; end
  phip = phipPrev + opts.relax*(phipSolved - phipPrev);

  if ~isempty(junctionList)
    [phin, phip] = junctionCouplingCorrection(dev, psi, phin, phip, V, Gprofile, junctionList);
  end

  outerStep = max([max(abs(psi-psiPrev)), max(abs(phin-phinPrev)), max(abs(phip-phipPrev))]);
  if opts.verbose
    [~,~,~,JnSegV,JpSegV] = assembleResidualQFPFD12(dev, psi, phin, phip, V, Gprofile);
    JtotV = JnSegV+JpSegV;
    fprintf('  outer=%d outerStep=%.4e Jmed=%.6e Jspread=%.4e\n', outerIter, outerStep, median(JtotV), max(JtotV)-min(JtotV));
  end
  if outerStep < opts.tolOuter && outerIter > 2
    [~, ~, ~, JnSegC, JpSegC] = assembleResidualQFPFD12(dev, psi, phin, phip, V, Gprofile);
    JtotC = JnSegC + JpSegC;
    JtermC = median(JtotC);
    deltaJabsC = max(JtotC) - min(JtotC);
    deltaJrelC = deltaJabsC / max(abs(JtermC), 1e-6);
    if deltaJabsC < opts.tolJabs && deltaJrelC < opts.tolJrel
      converged = true;
      break;
    end
  end
end

[F, n, p, JnSeg, JpSeg] = assembleResidualQFPFD12(dev, psi, phin, phip, V, Gprofile);
isFiniteState = isFiniteState && all(isfinite(F)) && all(isfinite(psi)) && ...
  all(isfinite(n)) && all(isfinite(p));

JtotSeg = JnSeg + JpSeg;
Jterminal = median(JtotSeg);
Jfloor = 1e-12;
deltaJ = (max(JtotSeg) - min(JtotSeg)) / max(max(abs(JtotSeg)), Jfloor);

idxP = 2:N-1;
chargeScale = c.Q * max([max(dev.Nd(idxP)), max(dev.Na(idxP)), max(n(idxP)), max(p(idxP)), 1e10]);
poissonResidual = max(abs(F(3*(idxP-1)+1))) / chargeScale;
poissonConverged = poissonResidual < 1e-6;

currentConservation = checkCurrentConservation(JtotSeg, Jterminal, 1e-7, 1e-3);
continuityResiduals = checkContinuityResiduals(dev, Gprofile, n, p, JnSeg, JpSeg);

chiRel = dev.chi - dev.chiRef;
EFn = c.VT*fd12inv(n./dev.Nc) - psi - chiRel;
EFp = -psi - chiRel - dev.Eg - c.VT*fd12inv(p./dev.Nv);
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
      equilibriumCheck.EFpSpan<1e-6 && qfpConsistency.pass;
  physicsPass = equilibriumCheck.pass;
else
  equilibriumCheck.pass = false;
  physicsPass = currentConservation.passAbs && currentConservation.passRel && ...
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

function [x, ok] = newtonSingleBlock(resFun, jacFun, x0, N, blockIdx, maxIter, tol, maxStep)
x = x0;
ok = true;
rowIdx = blockIdx:3:3*N;
for iter = 1:maxIter
  Ffull = resFun(x);
  Fblock = Ffull(rowIdx);
  if ~all(isfinite(Fblock))
    ok = false; return;
  end
  jacFunFull = @(xv) subset(jacFun(xv), rowIdx);
  J = assembleJacobianFDColored(jacFunFull, x, N, 1, 1);   % coloured FD Jacobian (identical to assembleJacobianFD.m, 6 instead of 2N evaluations)
  delta = J \ (-Fblock);
  if ~all(isfinite(delta))
    ok = false; return;
  end
  delta = max(min(delta, maxStep), -maxStep);
  normF0 = norm(Fblock, 2);
  alpha = 1.0;
  bestAlpha = 0; bestNorm = normF0;
  for bt = 1:8
    xTrial = x + alpha*delta;
    Ftrial = resFun(xTrial);
    Ftb = Ftrial(rowIdx);
    if all(isfinite(Ftb))
      nF = norm(Ftb, 2);
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
  x = x + step;
  if max(abs(step)) < tol && iter > 1
    break;
  end
end
end

function v = subset(Ffull, rowIdx)
v = Ffull(rowIdx);
end

function [phin, phip] = junctionCouplingCorrection(dev, psi, phin, phip, V, Gprofile, junctionList)
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
Ffull = assembleResidualQFPFD12(dev, psi, phinT, phipT, V, Gprofile);
F = Ffull(rows);
end
