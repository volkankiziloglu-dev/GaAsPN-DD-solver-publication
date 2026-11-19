function out = gummelSolveQFP(dev, V, Gprofile, opts)
%GUMMELSOLVEQFP Gummel-alternating outer loop over QFP variables.
%
% Combines the two things this project's Newton-only QFP solver
% (solveCoupledQFP.m) and the pre-existing Gummel/thomasSolve pipeline
% each got half of, matching how SCAPS-1D and AFORS-HET actually solve
% this class of problem (both use Gummel's alternating scheme with
% Newton substeps per equation block, NOT one fully-coupled system --
% confirmed against public descriptions of both tools):
%   - carrierDensitiesQFP.m's n,p = Nc*exp(...)/Nv*exp(...) are
%     structurally positive for any finite argument -- the actual N.16
%     fix, kept from solveCoupledQFP.m.
%   - Each Gummel sub-problem (Poisson in psi alone; electron continuity
%     in phin alone; hole continuity in phip alone) is solved with the
%     OTHER TWO variable blocks frozen -- so no single linear solve ever
%     has to mix a majority-carrier row/column scale with a
%     minority-carrier one in the same matrix. This is exactly the
%     structural problem diagnosed in solveCoupledQFP.m's fully-coupled
%     Newton (a Jacobian spanning ~80 orders of magnitude once JtJ was
%     formed for Levenberg-Marquardt; a raw Newton step of ~1e11 V in a
%     single poorly-determined direction) -- Gummel alternation sidesteps
%     it by construction, not by a smarter linear-algebra trick.
%
% Each sub-problem is itself nonlinear (n,p depend exponentially on the
% frozen-vs-solved variables) and is solved by its OWN small Newton
% iteration, reusing assembleResidualQFP.m/assembleJacobianFD.m with
% varsPerNode=1 -- the resulting per-block Jacobian is exactly
% tridiagonal (each row only couples to its immediate spatial
% neighbors IN THE SAME VARIABLE), well-scaled (O(1) volt-scale, since
% within one block there is no majority/minority disparity), and solved
% with MATLAB's sparse \ (no thomasSolve.m pivoting risk reintroduced).
c = constants();
N = dev.N;

if nargin < 4, opts = struct(); end
if ~isfield(opts,'maxOuter') || isempty(opts.maxOuter), opts.maxOuter = 80; end
if ~isfield(opts,'maxInner') || isempty(opts.maxInner), opts.maxInner = 30; end
if ~isfield(opts,'tolOuter') || isempty(opts.tolOuter), opts.tolOuter = 1e-10; end
if ~isfield(opts,'tolInner') || isempty(opts.tolInner), opts.tolInner = 1e-13; end
if ~isfield(opts,'relax') || isempty(opts.relax), opts.relax = 1.0; end
if ~isfield(opts,'verbose') || isempty(opts.verbose), opts.verbose = false; end
% opts.rowScale (default false, dd_solver v1.1.2): scale each block row by its largest Jacobian entry (Jacobian, residual and
% line-search norm alike). Off by default, so existing results are unchanged.
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

% N.16.18: junction node pairs needing the local fully-coupled correction
% (see junctionCouplingCorrection below) -- built once, reused every outer
% iteration. Empty when interface recombination is off (Sn_int/Sp_int<=0
% or dev.Sn_int_j/Sp_int_j all zero), matching interfaceRecombination.m's
% own "returns zero, ideal interface" fallback exactly, so devices without
% interface recombination take an identical path to before this fix.
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

  % Gprofile is dropped to 0 ONLY for the Jacobian probe passed to
  % newtonSingleBlock, never for the residual F itself -- same fix as
  % solveCoupledQFP.m: Gprofile enters assembleResidualQFP purely
  % additively (never multiplied by psi/phin/phip), so dF/dx is
  % mathematically identical whether G is the real illumination profile
  % or zero. Left in during FD probing, a large G (~1e21-1e22 under
  % AM1.5G) swamps the double-precision representation of F and the
  % derivative of the much smaller recombination term comes out exactly
  % zero -- reproduced here in the Gummel-alternating solver too
  % (isFinite=0 at the very first illuminated point before this fix).

  % --- Sub-problem 1: Poisson in psi alone (phin,phip frozen) ---
  resFunPsi = @(psiVec) assembleResidualQFP(dev, psiVec, phin, phip, V, Gprofile);
  jacFunPsi = @(psiVec) assembleResidualQFP(dev, psiVec, phin, phip, V, Gzero);
  [psiSolved, okPsi] = newtonSingleBlock(resFunPsi, jacFunPsi, psi, N, 1, opts.maxInner, opts.tolInner, 0.5, opts.rowScale);
  if ~okPsi, isFiniteState = false; break; end
  % Classic Gummel under-relaxation: blend toward the sub-solved value
  % rather than accepting it outright. Diagnosed directly on this device
  % (illuminated, V=-0.05): full-step (opts.relax=1) Gummel converged the
  % right ORDER OF MAGNITUDE (J~1.2e-2 vs target 1.38e-2 A/cm^2, isFinite)
  % but then sat at an EXACTLY unchanged fixed point for 150 outer
  % iterations with current still not spatially conserved
  % (deltaJabs~1.3e-2) -- a textbook undamped-Gummel oscillation/stall,
  % not a slow-convergence issue (more iterations made zero difference).
  psi = psiPrev + opts.relax*(psiSolved - psiPrev);
  psi(1) = dev.psiEq(1);
  psi(N) = dev.psiEq(N) + V;

  % --- Sub-problem 2: electron continuity in phin alone (psi,phip frozen) ---
  resFunPhin = @(phinVec) assembleResidualQFP(dev, psi, phinVec, phip, V, Gprofile);
  jacFunPhin = @(phinVec) assembleResidualQFP(dev, psi, phinVec, phip, V, Gzero);
  [phinSolved, okPhin] = newtonSingleBlock(resFunPhin, jacFunPhin, phin, N, 2, opts.maxInner, opts.tolInner, 1.0, opts.rowScale);
  if ~okPhin, isFiniteState = false; break; end
  phin = phinPrev + opts.relax*(phinSolved - phinPrev);

  % --- Sub-problem 3: hole continuity in phip alone (psi,phin frozen) ---
  resFunPhip = @(phipVec) assembleResidualQFP(dev, psi, phin, phipVec, V, Gprofile);
  jacFunPhip = @(phipVec) assembleResidualQFP(dev, psi, phin, phipVec, V, Gzero);
  [phipSolved, okPhip] = newtonSingleBlock(resFunPhip, jacFunPhip, phip, N, 3, opts.maxInner, opts.tolInner, 1.0, opts.rowScale);
  if ~okPhip, isFiniteState = false; break; end
  phip = phipPrev + opts.relax*(phipSolved - phipPrev);

  % --- Sub-problem 4 (N.16.18): local fully-coupled correction at each
  % interface-recombination junction. interfaceRecombination.m couples
  % n(iE),n(iA),p(iE),p(iA) (both carriers, both nodes straddling the
  % junction) into ONE recombination expression feeding all four of
  % Fphin(iE)/Fphin(iA)/Fphip(iE)/Fphip(iA) -- a correction direction the
  % three single-variable blocks above cannot express (each freezes the
  % other two variables, so phin's block can't see phip(iE)/phip(iA)
  % moving together with it, and vice versa). Confirmed by direct
  % instrumentation (N.16.16): without this, massActionMaxRel is an EXACT
  % fixed point independent of outer iteration count (2000 iterations
  % changed nothing to the last digit) for any nonzero interface
  % recombination velocity, however small. This local 4-unknown Newton
  % step (dense FD Jacobian, cheap) restores exactly the coupling the
  % block decomposition is structurally missing, without reintroducing
  % solveCoupledQFP.m's device-wide conditioning catastrophe (only 2
  % nodes x 2 carriers here, not 2 nodes x 2 carriers x N).
  if ~isempty(junctionList)
    [phin, phip] = junctionCouplingCorrection(dev, psi, phin, phip, V, Gprofile, junctionList);
  end

  outerStep = max([max(abs(psi-psiPrev)), max(abs(phin-phinPrev)), max(abs(phip-phipPrev))]);
  if opts.verbose
    [~,nV,pV,JnSegV,JpSegV] = assembleResidualQFP(dev, psi, phin, phip, V, Gprofile);
    JtotV = ifaceCorrectTotalCurrent(dev, JnSegV+JpSegV, nV, pV);
    fprintf('  outer=%d outerStep=%.4e Jmed=%.6e Jspread=%.4e\n', outerIter, outerStep, median(JtotV), max(JtotV)-min(JtotV));
  end
  % N.16.13 false-convergence fix: step-size alone is not sufficient to
  % declare the outer Gummel loop converged (a stalled fixed point has
  % outerStep=0 too, see solveElectronContinuity comment above). Also
  % require the SAME combined abs+rel current-conservation test the
  % acceptance test applies downstream, checked directly on this
  % iteration's Jtot(x) -- if it doesn't hold, keep iterating (or exhaust
  % maxOuter and report converged=false, which is now the honest answer).
  if outerStep < opts.tolOuter && outerIter > 2
    [~, nC, pC, JnSegC, JpSegC] = assembleResidualQFP(dev, psi, phin, phip, V, Gprofile);
    JtotC = ifaceCorrectTotalCurrent(dev, JnSegC + JpSegC, nC, pC);   % unchanged unless 'four'/'cross' (dd_solver v1.1.1)
    JtermC = median(JtotC);
    deltaJabsC = max(JtotC) - min(JtotC);
    deltaJrelC = deltaJabsC / max(abs(JtermC), 1e-6);
    if deltaJabsC < opts.tolJabs && deltaJrelC < opts.tolJrel
      converged = true;
      break;
    end
  end
end

[F, n, p, JnSeg, JpSeg] = assembleResidualQFP(dev, psi, phin, phip, V, Gprofile);
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
continuityResiduals = checkContinuityResiduals(dev, Gprofile, n, p, JnSeg, JpSeg);

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
  % N.16.13: strict AND, not the permissive OR this used to use (which let
  % e.g. a passAbs-only or continuityResiduals-only pass through). Matches
  % the acceptance test's own combined-criterion definition exactly, so
  % out.physicsPass is trustworthy on its own without re-deriving it.
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

function [x, ok] = newtonSingleBlock(resFun, jacFun, x0, N, blockIdx, maxIter, tol, maxStep, rowScale)
% Newton solve for ONE variable block (blockIdx in {1=psi,2=phin,3=phip}
% of assembleResidualQFP's 3-per-node packing), all others held fixed by
% resFun's/jacFun's own closure. The corresponding row range of F is
% exactly the rows for that block; the Jacobian (varsPerNode=1 call to
% assembleJacobianFD) is exactly tridiagonal since each row only depends
% on its immediate spatial neighbors WITHIN the same variable. jacFun is
% the SAME residual but with Gprofile zeroed (see caller's comment) --
% used only to build the Jacobian, never to evaluate the actual F/merit.
%
% N.16.13: an earlier version of this fix made `ok` FAIL immediately
% whenever a single call's line search stagnated (bestAlpha=0) or
% exhausted maxIter without reaching a strict residual tolerance, on the
% reasoning that step-size-only convergence (the ORIGINAL logic here)
% could mistake "stuck at a large residual" for "converged". Verified by
% direct instrumentation this DOES matter, but gating individual calls
% on it is the wrong layer to fix it at: at V=0.3-1.2 (previously-clean
% PASSES) the phin/phip blocks routinely make partial progress in ONE
% outer iteration's call (e.g. residual 1.93e-3 -> 1.83e-3) and then
% genuinely stagnate for THIS call -- not because the overall problem is
% unconverged, but because Gummel alternation is block-coordinate by
% design: the NEXT outer iteration's call sees updated psi/phip and
% typically continues making progress from there. Hard-failing on a
% single call's stagnation aborted the outer loop after 1-2 iterations
% even for these known-good points. The false-convergence problem is
% real but belongs at the OUTER level (see gummelSolveQFP's outerStep/
% current-conservation combined check, and physicsPass's strict-AND
% below) where the actual multi-pass residual is what's evaluated --
% `ok` here stays a pure numerical-sanity (isfinite) gate, matching what
% it always effectively was, so the outer loop's own iteration
% mechanism is not disrupted by this fix.
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
  F0forJ = jacFunFull(x);
  J = assembleJacobianFD(jacFunFull, x, F0forJ, N, 1, 1);
  if rowScale
    weights=1./max(full(max(abs(J),[],2)),realmin);
    delta=(spdiags(weights,0,N,N)*J)\(-weights.*Fblock);
  else
    weights=ones(N,1);delta=J\(-Fblock);
  end
  if ~all(isfinite(delta))
    ok = false; return;
  end
  delta = max(min(delta, maxStep), -maxStep);
  normF0 = norm(weights.*Fblock, 2);
  alpha = 1.0;
  bestAlpha = 0; bestNorm = normF0;
  for bt = 1:8
    xTrial = x + alpha*delta;
    Ftrial = resFun(xTrial);
    Ftb = Ftrial(rowIdx);
    if all(isfinite(Ftb))
      nF = norm(weights.*Ftb, 2);
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
% N.16.18: local fully-coupled Newton correction for the 4 unknowns
% {phin(iE),phin(iA),phip(iE),phip(iA)} at each interface-recombination
% junction (iE,iA straddle the junction -- see interfaceRecombination.m).
% Dense 4x4 FD Jacobian, small damped Newton loop; psi and every OTHER
% node's phin/phip stay fixed (this is a LOCAL correction, not a second
% device-wide coupled solve). No-op (returns unchanged) if this junction's
% residual is already tiny, so it costs nothing once converged.
for jk = 1:size(junctionList,1)
  iE = junctionList(jk,1); iA = junctionList(jk,2);
  rows = [3*(iE-1)+2, 3*(iA-1)+2, 3*(iE-1)+3, 3*(iA-1)+3]; % Fphin(iE),Fphin(iA),Fphip(iE),Fphip(iA)
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
    delta = max(min(delta, 0.5), -0.5); % damp, matches phin/phip block maxStep
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
Ffull = assembleResidualQFP(dev, psi, phinT, phipT, V, Gprofile);
F = Ffull(rows);
end
