function out = gummelSolveML(dev, V, Gprofile, opts)
% Copy of RearHJT/gummelSolveML.m -- calls newtonPoissonML (eps0-corrected).
% The electron/hole BC dispatch was already conditional (RearHJT's own
% "Item 6" fix, missing in common/gummelSolve.m) -- unchanged
% here. Default relax/maxIter were changed with the SAME logic as in gummelSolve.m
% (see the note in that file).
%
% DIAGNOSTIC ADDITIONS (external review items 1.2/1.5/3.3):
%  - the exp() arguments in the un/vp computation are now clamped
%    (the same kind of protection sgBernoulli already applies) -- with large band
%    offsets/forward bias there was a risk of overflow/NaN.
%  - newtonPoissonML now returns [psi,converged,iter,residual];
%    this information is passed out as out.poissonConverged/out.poissonResidual,
%    and the outer Gummel loop now ALSO includes a failed inner Poisson solve
%    in the "converged" status (previously it was silently ignored).
%  - the current is computed NOT ONLY at floor(N/2) but at ALL interior nodes, and its
%    continuity (deltaJ) is reported.
c = constants();
N = dev.N; ni = dev.ni; Nd = dev.Nd; Na = dev.Na;

if nargin < 4, opts = struct(); end
if ~isfield(opts,'maxIter'), opts.maxIter = 1500; end
if ~isfield(opts,'relax'), opts.relax = 0.05; end

EXP_CLAMP = 700; % exp(700)~1e304, safe limit before double overflow (~1.8e308)

if Nd(1) > 0
  nLeft = Nd(1);
else
  nLeft = ni(1)^2 / Na(1);
end
pLeft = ni(1)^2 / nLeft;
if Na(N) > 0
  pRight = Na(N);
else
  pRight = ni(N)^2 / Nd(N);
end
nRight = ni(N)^2 / pRight;

if isfield(opts,'psiInit') && ~isempty(opts.psiInit)
  psi = opts.psiInit;
  psiRightTarget = dev.psiEq(N) + V;
else
  psi = dev.psiEq;
  psiRightTarget = dev.psiEq(N) + V;
  psi(N) = psiRightTarget;
end
if isfield(opts,'nInit') && ~isempty(opts.nInit), n = opts.nInit; else, n = dev.nEq; end
if isfield(opts,'pInit') && ~isempty(opts.pInit), p = opts.pInit; else, p = dev.pEq; end

Sback = dev.Sback;
converged = false;
poissonConverged = true; % status of the last inner Poisson call
poissonResidual = NaN;
poissonFailCount = 0;
negFloorTotal = 0;
negFloorMaxPerIter = 0;
subFloorPositiveTotal = 0;
minNRaw = Inf;
minPRaw = Inf;
for iter = 1:opts.maxIter
  nPrevIter = n; pPrevIter = p;

  if Nd(N) > 0
    nNew = solveElectronContinuity(dev, psi, nLeft, nRight, Gprofile, nPrevIter, pPrevIter);
  else
    nNew = solveElectronContinuity(dev, psi, nLeft, nRight, Gprofile, nPrevIter, pPrevIter, Sback, nRight);
  end
  if Na(N) > 0
    pNew = solveHoleContinuity(dev, psi, pLeft, pRight, Gprofile, nPrevIter, pPrevIter);
  else
    pNew = solveHoleContinuity(dev, psi, pLeft, pRight, Gprofile, nPrevIter, pPrevIter, Sback, pRight);
  end
  % N.16 (external review): the 1e-30 floor previously SILENTLY made negative/meaningless
  % solutions (destructive-cancellation failure of thomasSolve.m, see
  % MATHEMATICAL_MODEL.md N.16) look "converged". The number of negative
  % (physically impossible, not merely positive-but-small) values
  % is counted separately and passed out as out.negFloorHit; this
  % count being nonzero at the last iteration makes out.converged FALSE
  % (below).
  negFloorCount = sum(nNew < 0) + sum(pNew < 0);
  negFloorTotal = negFloorTotal + negFloorCount;
  negFloorMaxPerIter = max(negFloorMaxPerIter, negFloorCount);
  subFloorPositiveTotal = subFloorPositiveTotal + ...
      sum(nNew >= 0 & nNew < 1e-30) + sum(pNew >= 0 & pNew < 1e-30);
  minNRaw = min(minNRaw, min(nNew));
  minPRaw = min(minPRaw, min(pNew));
  nNew(nNew < 1e-30) = 1e-30;
  pNew(pNew < 1e-30) = 1e-30;

  maxRelNP = max(max(abs(log(nNew./n))), max(abs(log(pNew./p))));
  n = nNew; p = pNew;

  psiEffN = psi + dev.psiOffsetN;
  psiEffP = psi + dev.psiOffsetP;
  argN = max(min(-psiEffN/c.VT, EXP_CLAMP), -EXP_CLAMP);
  argP = max(min(psiEffP/c.VT, EXP_CLAMP), -EXP_CLAMP);
  un = n .* exp(argN);
  vp = p .* exp(argP);
  [psiNew, poissonConverged, ~, poissonResidual] = newtonPoissonML(dev, psi, un, vp, dev.psiEq(1), psiRightTarget, false);
  if ~poissonConverged
    poissonFailCount = poissonFailCount + 1;
  end

  step = opts.relax * (psiNew - psi);
  psi = psi + step;
  maxChange = max(abs(step));

  if maxChange < 1e-12 && maxRelNP < 1e-9 && iter > 3 && poissonConverged
    converged = true;
    break;
  end
end

% The loop updates psi after solving n/p. At low illumination a single
% carrier re-solve is not enough because SRH recombination couples the two
% continuity equations. Iterate them at fixed final psi until both carrier
% profiles are mutually consistent.
negFloorCount = 0;
for polishIter = 1:25
  nBeforePolish = n;
  pBeforePolish = p;
  if Nd(N) > 0
    nPolish = solveElectronContinuity(dev, psi, nLeft, nRight, Gprofile, n, p);
  else
    nPolish = solveElectronContinuity(dev, psi, nLeft, nRight, Gprofile, n, p, Sback, nRight);
  end
  nPolishRaw = nPolish;
  nPolish = max(nPolishRaw, 1e-30);
  if Na(N) > 0
    pPolish = solveHoleContinuity(dev, psi, pLeft, pRight, Gprofile, nPolish, p);
  else
    pPolish = solveHoleContinuity(dev, psi, pLeft, pRight, Gprofile, nPolish, p, Sback, pRight);
  end
  negPolish = sum(nPolishRaw < 0) + sum(pPolish < 0);
  negFloorTotal = negFloorTotal + negPolish;
  negFloorCount = negFloorCount + negPolish;
  negFloorMaxPerIter = max(negFloorMaxPerIter, negPolish);
  subFloorPositiveTotal = subFloorPositiveTotal + ...
      sum(nPolishRaw >= 0 & nPolishRaw < 1e-30) + sum(pPolish >= 0 & pPolish < 1e-30);
  minNRaw = min(minNRaw, min(nPolishRaw));
  minPRaw = min(minPRaw, min(pPolish));
  n = max(nPolish, 1e-30);
  p = max(pPolish, 1e-30);
  polishChange = max(max(abs(log(n./nBeforePolish))), ...
                     max(abs(log(p./pBeforePolish))));
  if polishChange < 1e-10
    break;
  end
end

% Current continuity: compute Jn+Jp at ALL interior nodes (not only floor(N/2)
% -- external review item 1.5). It should be constant in steady state;
% if not, deltaJ captures it.
idxJ = 1:N-1;
psiEffFinalN = psi + dev.psiOffsetN;
psiEffFinalP = psi + dev.psiOffsetP;
dPsiN = (psiEffFinalN(idxJ+1)-psiEffFinalN(idxJ))/c.VT;
dPsiP = (psiEffFinalP(idxJ+1)-psiEffFinalP(idxJ))/c.VT;
DnSeg = c.VT*0.5*(dev.muN(idxJ)+dev.muN(idxJ+1));
DpSeg = c.VT*0.5*(dev.muP(idxJ)+dev.muP(idxJ+1));
JnSeg = c.Q*DnSeg./dev.h(idxJ) .* (n(idxJ+1).*sgBernoulli(dPsiN) - n(idxJ).*sgBernoulli(-dPsiN));
JpSeg = c.Q*DpSeg./dev.h(idxJ) .* (p(idxJ).*sgBernoulli(dPsiP) - p(idxJ+1).*sgBernoulli(-dPsiP));
JtotSeg = JnSeg + JpSeg;

Jterminal = median(JtotSeg);
Jfloor = 1e-12; % A/cm^2, absolute floor (prevents division of deltaJ by zero)
deltaJ = (max(JtotSeg) - min(JtotSeg)) / max(max(abs(JtotSeg)), Jfloor);

isFinite = all(isfinite(psi)) && all(isfinite(n)) && all(isfinite(p)) && isfinite(Jterminal);
% At indoor injection levels the raw-density formulation reaches its
% floating-point conditioning floor near the highly doped front contact.
% Keep a strict relative test, but permit at most 0.1 uA/cm^2 absolute
% imbalance; this remains negligible for outdoor currents and is reported.
currentConservation = checkCurrentConservation(JtotSeg, Jterminal, 1e-7, 1e-3);
continuityResiduals = checkContinuityResiduals(dev, Gprofile, n, p, JnSeg, JpSeg);
chiRel = dev.chi-dev.chiRef;
EFn = c.VT*log(n./dev.Nc)-psi-chiRel;
EFp = -psi-chiRel-dev.Eg-c.VT*log(p./dev.Nv);
equilibriumCheck = struct();
equilibriumCheck.EFnSpan = max(EFn)-min(EFn);
equilibriumCheck.EFpSpan = max(EFp)-min(EFp);
equilibriumCheck.massActionMaxRel = max(abs(n.*p-dev.ni.^2)./max(dev.ni.^2,1e-300));
equilibriumLike = abs(V)<1e-14 && max(abs(Gprofile))<1;
if equilibriumLike
  equilibriumCheck.pass = equilibriumCheck.EFnSpan<1e-6 && ...
      equilibriumCheck.EFpSpan<1e-6 && equilibriumCheck.massActionMaxRel<1e-6;
  physicsPass = negFloorTotal==0 && currentConservation.pass && equilibriumCheck.pass;
else
  equilibriumCheck.pass = false;
  physicsPass = negFloorTotal==0 && currentConservation.pass && ...
      (continuityResiduals.pass || currentConservation.passAbs);
end

out.psi = psi; out.n = n; out.p = p; out.J = Jterminal;
% N.16 (external review, IMPORTANT FIX): in a first attempt, setting
% out.converged=false when negFloorHit>0 was tried AND REVERTED -- it broke the illuminated N.1
% FULL sweepJVML.m sweep immediately, EVEN at V=0.0 (see
% MATHEMATICAL_MODEL.md N.16.8) -- this revealed that the defect diagnosed by N.16
% is MUCH more widespread than assumed (not specific to high V; some
% negative-floored nodes appear at almost EVERY processed point).
% Instead of silently putting a behavior change of this magnitude
% into production, negFloorHit is reported FOR DIAGNOSIS ONLY --
% out.converged stays UNCHANGED with its current definition (isFinite+convergence),
% and calling code (such as sweepJVML.m) can use this information IF IT WANTS.
out.negFloorHit = negFloorCount;
out.negFloorTotal = negFloorTotal;
out.negFloorMaxPerIter = negFloorMaxPerIter;
out.subFloorPositiveTotal = subFloorPositiveTotal;
out.minNRaw = minNRaw;
out.minPRaw = minPRaw;
out.currentConservation = currentConservation;
out.continuityResiduals = continuityResiduals;
out.equilibriumCheck = equilibriumCheck;
out.physicsPass = physicsPass;
out.converged = converged && isFinite && physicsPass;
out.iters = iter;
out.isFinite = isFinite;
out.poissonConverged = poissonConverged;
out.poissonResidual = poissonResidual;
out.poissonFailCount = poissonFailCount;
out.deltaJ = deltaJ;
out.Jseg = JtotSeg;
end
