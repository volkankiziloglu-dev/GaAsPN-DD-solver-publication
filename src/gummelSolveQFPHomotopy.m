function out = gummelSolveQFPHomotopy(dev, V, Gprofile, opts)
%GUMMELSOLVEQFPHOMOTOPY Adaptive illumination-ramp continuation around
% gummelSolveQFP.m, targeting the V~0 phi_n-stall (N.16.10/N.16.13/N.16.14)
% at its actual mechanism: jumping straight from a dark/low-injection
% initial guess to full illumination in one Gummel solve asks the
% electron continuity block to absorb a huge injection-profile change in
% one step. Instead, ramp lambda in G_lambda = lambda*Gprofile from 0 to
% 1, warm-starting each step from the previous ACCEPTED step's
% (psi,phin,phip) -- halving the step on failure, growing it on repeated
% success, giving up only once the step is below a floor.
%
% Step-acceptance during the ramp uses isFinite && physicsPass (the
% direct physics-correctness check: current conservation + continuity
% residuals), not the full qfpResultAccepted.m gate -- physicsPass is
% what actually matters for "is this a legitimate DC operating point to
% warm-start the next lambda from", while qfpResultAccepted's additional
% r.converged requirement (outer step size below tolOuter) is a stricter
% numerical-tidiness criterion that need not hold at every intermediate
% lambda. The FINAL returned result's out.accepted field uses the full
% qfpResultAccepted gate, matching every other script in this project.
N = dev.N;

if nargin < 4, opts = struct(); end
if ~isfield(opts,'dLambdaInit') || isempty(opts.dLambdaInit), opts.dLambdaInit = 0.1; end
if ~isfield(opts,'dLambdaMax') || isempty(opts.dLambdaMax), opts.dLambdaMax = 0.25; end
if ~isfield(opts,'dLambdaMin') || isempty(opts.dLambdaMin), opts.dLambdaMin = 1e-4; end
if ~isfield(opts,'growFactor') || isempty(opts.growFactor), opts.growFactor = 1.5; end
if ~isfield(opts,'verbose') || isempty(opts.verbose), opts.verbose = false; end
if ~isfield(opts,'psiInit') || isempty(opts.psiInit), opts.psiInit = dev.psiEq; end
if ~isfield(opts,'phinInit') || isempty(opts.phinInit), opts.phinInit = zeros(N,1); end
if ~isfield(opts,'phipInit') || isempty(opts.phipInit), opts.phipInit = zeros(N,1); end

% Log-scale step schedule (opt-in, default OFF -- the existing linear
% dLambda schedule above is untouched and remains the default for every
% caller that doesn't ask for this). Targets the indoor (dim-light)
% convergence regime: a benchmark under 2700K/1000lux found the LINEAR
% schedule's first step (dark -> lambda=dLambdaInit=0.1, a full 10% of
% full illumination in one jump) never converges even after >29 minutes,
% because indoor's absolute generation rate is ~1e3-1e4x smaller than
% AM1.5G's, so a "10%" step is still a huge RELATIVE jump in the
% low-injection regime the electron continuity block is most sensitive
% to. Geometric-in-lambda stepping (lambda multiplied, not added, each
% successful step) naturally produces tiny absolute steps near lambda=0
% and only grows once further from it -- matching where the difficulty
% actually lives. See MATHEMATICAL_MODEL.md N.16.21/N.17.
if ~isfield(opts,'useLogSchedule') || isempty(opts.useLogSchedule), opts.useLogSchedule = false; end
if ~isfield(opts,'lambdaLogStart') || isempty(opts.lambdaLogStart), opts.lambdaLogStart = 1e-8; end
if ~isfield(opts,'dLogLambdaInit') || isempty(opts.dLogLambdaInit), opts.dLogLambdaInit = 0.5; end
if ~isfield(opts,'dLogLambdaMax') || isempty(opts.dLogLambdaMax), opts.dLogLambdaMax = 1.0; end
if ~isfield(opts,'dLogLambdaMin') || isempty(opts.dLogLambdaMin), opts.dLogLambdaMin = 1e-3; end

% Inner per-lambda solver (opt-in, default @gummelSolveQFP -- unchanged
% behavior for every existing caller). gummelSolveQFPAnalytic.m is bit-
% identical to gummelSolveQFP.m (cross-validated Jacobians, see
% assemblePsi/Phin/PhipJacobianAnalyticQFP.m and
% tests/test_QFP_analytic_jacobian_crosscheck.m) but ~100x faster per
% call on realistic device sizes -- since the homotopy ramp's cost is
% dominated by repeated per-lambda solves (each a full outer Gummel
% loop, each up to opts internal maxOuter), passing
% @gummelSolveQFPAnalytic here directly attacks the OTHER half of
% indoor's slowness (per-call cost), independent of the log-schedule fix
% above (which attacks wasted/failed steps). The two are orthogonal and
% meant to be used together for indoor.
if ~isfield(opts,'innerSolver') || isempty(opts.innerSolver), opts.innerSolver = @gummelSolveQFP; end
innerSolver = opts.innerSolver;

lambda = 0;
dLambda = opts.dLambdaInit;
dLogLambda = opts.dLogLambdaInit;
% logPos: current position in log10(lambda) space, tracked SEPARATELY
% from lambda so the very-first trial (lambda still 0, log10(0)=-Inf
% undefined) can shrink on failure exactly like every later step -- an
% earlier version of this schedule re-tried the SAME fixed
% opts.lambdaLogStart on every failure at lambda=0 (since it read
% lambda<=0 as a special case rather than going through logPos), which
% for a device where even that first tiny lambda is too big a jump just
% burned through the failure floor doing nothing. Seeded so the first
% trial (logPos+dLogLambda) equals log10(lambdaLogStart) exactly.
logPos = log10(opts.lambdaLogStart) - opts.dLogLambdaInit;
psi = opts.psiInit; phin = opts.phinInit; phip = opts.phipInit;

% lambda=0 anchor: dark (or, if V!=0, unilluminated bias) solve. Gated on
% isFinite ONLY, not physicsPass -- physicsPass's equilibriumLike branch
% (V==0 && G==0) demands massActionMaxRel<1e-6, a MUCH tighter bar than
% what the ramp itself needs. Found by direct instrumentation: MD421's
% own dev.psiEq (from the separate solveEquilibrium.m, not gummelSolveQFP)
% is not perfectly self-consistent for this device's more complex
% structure (massActionMaxRel=0.238 observed, i.e. 24% -- solveEquilibrium.m
% itself under-converges for this specific 3-layer/2-junction/multi-center-
% SRH device), yet the SAME starting point is a perfectly good warm-start
% seed for the illuminated ramp -- gating the anchor on physicsPass made
% the homotopy refuse to even begin for MD421 (lambdaReached stuck at 0
% for every monochromatic energy tested), even though the direct
% (non-homotopy) illuminated V=0 solve from the SAME dev.psiEq converges
% fine via the outer Gummel loop's own iteration. Newton-based homotopy
% exists precisely to improve on an imperfect seed; refusing to start
% from one defeats the purpose.
r = innerSolver(dev, V, 0*Gprofile, struct('psiInit',psi,'phinInit',phin,'phipInit',phip));
history = struct('lambda',{0},'accepted',{r.isFinite},'dLambda',{dLambda});
if ~r.isFinite
  out = r;
  out.homotopy = struct('lambdaReached',0,'homotopyConverged',false,'history',history);
  out.accepted = false;
  return;
end
psi = r.psi; phin = r.phin; phip = r.phip;

while lambda < 1 - 1e-12
  if opts.useLogSchedule
    lambdaTrial = min(1, 10^(logPos + dLogLambda));
  else
    lambdaTrial = min(1, lambda + dLambda);
  end
  rTrial = innerSolver(dev, V, lambdaTrial*Gprofile, struct('psiInit',psi,'phinInit',phin,'phipInit',phip));
  stepOk = rTrial.isFinite && rTrial.physicsPass;
  if opts.verbose
    fprintf('  homotopy lambda %.4e -> %.4e (dLambda=%.4e dLogLambda=%.4e) stepOk=%d\n', lambda, lambdaTrial, dLambda, dLogLambda, stepOk);
  end
  history(end+1) = struct('lambda',lambdaTrial,'accepted',stepOk,'dLambda',dLambda); %#ok<AGROW>
  if stepOk
    lambda = lambdaTrial;
    if opts.useLogSchedule && lambda > 0
      logPos = log10(lambda); % anchor to the actual accepted value, avoid drift
    end
    psi = rTrial.psi; phin = rTrial.phin; phip = rTrial.phip;
    r = rTrial;
    dLambda = min(opts.dLambdaMax, dLambda*opts.growFactor);
    dLogLambda = min(opts.dLogLambdaMax, dLogLambda*opts.growFactor);
  else
    dLambda = dLambda/2;
    dLogLambda = dLogLambda/2;
    if opts.useLogSchedule
      if dLogLambda < opts.dLogLambdaMin
        break;
      end
    else
      if dLambda < opts.dLambdaMin
        break;
      end
    end
  end
end

out = r;
out.homotopy = struct('lambdaReached',lambda,'homotopyConverged',lambda>=1-1e-9,'history',history);
out.accepted = qfpResultAccepted(r) && out.homotopy.homotopyConverged;
end
