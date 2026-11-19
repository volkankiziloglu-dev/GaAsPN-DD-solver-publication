function [pts, status] = sweepJVML(dev, Gprofile, Vmin, Vmax, Vstep, initialState)
% Copy of RearHJT/sweepJVML.m -- calls gummelSolveML (eps0-corrected),
% PLUS a V=0->Vmin warm-up sweep (after the eps0 fix, if
% Vmin is far from zero, a single large Newton jump could lead to catastrophic
% cancellation/NaN in thomasSolve.m -- see MATHEMATICAL_MODEL.md
% Section K items 12/13).
%
% Monotonic-J stopping threshold is RELATIVE (fix of the absolute-tolerance
% error found in N.7): JtolAbs = max(1e-7, 1e-3*J0).
%
% FIX (external review items 1.1/1.4): r.converged && r.isFinite is checked at EVERY
% point; on failure, starting from the previous converged
% point, the voltage step is halved and retried (continuation/
% bisection, at most maxSubdiv levels).
%
% FIX (external review item 2.4, REVISED): the originally proposed
% idea "add deltaJ to the acceptance threshold" was TRIED and shown to be
% WRONG -- with debug_deltaJ_calibration.m/debug_deltaJ_profile_good.m
% it was found that EVEN the ALREADY validated N.1 homojunction under illumination
% gives deltaJ~1.4-1.8 (%140-180!) -- LARGER than the "bad"
% ~0.1-0.4 of Config A. Reason: under illumination Jtot(x) is not expected
% to be FLAT physically in the first place (the Beer-Lambert absorption profile changes G(x)
% with depth, so the total current at each point varies with the local
% generation-recombination budget -- it is expected to be flat ONLY in the
% dark). In this form deltaJ cannot distinguish good from bad.
% Instead, the quantity that REALLY catches Config A is used:
% the photocurrent PHYSICALLY CANNOT EXCEED the MAXIMUM absorbed photon flux (JscMax)
% -- a solid conservation bound, not an arbitrary threshold. This bound
% is applied ONLY in the J>0 (photocurrent-dominated, near Jsc) regime; once V
% exceeds Voc, J becomes strongly negative (dark-current-dominated) and
% |J| may exceed JscMax -- this is NORMAL and is not limited.
%
% A full generation/recombination-budget-aware continuity-residual
% check (what the external review actually asked for) was deferred as a SEPARATE, larger
% task.
%
% FIX (external review item 2.6): a `status` struct was added as a second (optional)
% output -- calling code can now see explicitly WHY the sweep
% stopped (normal Vmax reached, no convergence, current exceeded the physical
% bound, or monotonic-J violation). Existing
% `pts = sweepJVML(...)` style old calls continue to work
% UNCHANGED (an extra output argument is backward compatible in MATLAB).
c = constants();
JscMax = c.Q * trapz(dev.x, Gprofile); % A/cm^2, physical upper bound in the J>0 regime
JscMaxMargin = 1.05; % %5 margin -- for coarse sweep/interpolation noise

status = struct('completed',false,'stopReason','','lastVoltage',NaN,'nPoints',0);

warmupStep = 0.002;
if nargin < 6 || isempty(initialState)
  prev = struct();
else
  required={'psi','n','p','J'};
  if ~all(isfield(initialState,required))
    error('sweepJVML:InitialStateFields', ...
      'initialState must contain psi, n, p, and J.');
  end
  prev=initialState;
end
prevV = 0;
if abs(Vmin) > warmupStep && nargin < 6
  Vwarm = 0:(warmupStep*sign(Vmin)):Vmin;
  if Vwarm(end) ~= Vmin, Vwarm(end+1) = Vmin; end
  for k = 1:numel(Vwarm)
    V = Vwarm(k);
    if isfield(prev,'psi')
      opts = struct('psiInit',prev.psi,'nInit',prev.n,'pInit',prev.p);
    else
      opts = struct();
    end
    r = gummelSolveML(dev, V, Gprofile, opts);
    if ~r.isFinite
      error('sweepJVML: NaN/Inf at V=%.4f in the warm-up sweep', V);
    end
    if ~r.converged
      error('sweepJVML: no convergence at V=%.4f in the warm-up sweep (iters=%d, poissonConverged=%d)', V, r.iters, r.poissonConverged);
    end
    if abs(V) < 1e-12 && abs(r.J) > JscMax*JscMaxMargin
      error('sweepJVML: in the warm-up sweep at V=%.4f, J=%.4e exceeds the physical JscMax=%.4e', V, r.J, JscMax);
    end
    prev = r;
    prevV = V;
  end
end
if ~isfield(prev,'psi')
  r = gummelSolveML(dev, Vmin, Gprofile, struct());
  if ~r.isFinite || ~r.converged
    error('sweepJVML: initial solution at V=%.4f did not converge/NaN (iters=%d)', Vmin, r.iters);
  end
  prev = r; prevV = Vmin;
end

Vs = Vmin:Vstep:Vmax;
n = numel(Vs);
pts = zeros(n,2);
nValid = 0;
Jprev = Inf;
J0 = abs(prev.J);
JtolAbs = max(1e-7, 1e-3*J0);
maxSubdiv = 4;

for k = 1:n
  Vtarget = Vs(k);
  if abs(Vtarget - prevV) < 1e-12
    r = prev; ok = true; failReason = '';
  else
    [r, ok, failReason] = stepToConverged(dev, Gprofile, prev, prevV, Vtarget, maxSubdiv, JscMax, JscMaxMargin);
  end
  if ~ok
    status.stopReason = failReason;
    status.lastVoltage = prevV;
    break;
  end
  if k > 1 && r.J > Jprev + JtolAbs
    status.stopReason = 'nonmonotonic';
    status.lastVoltage = prevV;
    break;
  end
  nValid = k;
  pts(k,:) = [Vtarget, r.J];
  Jprev = r.J;
  prev = r; prevV = Vtarget;
end
pts = pts(1:nValid,:);
status.nPoints = nValid;
if isempty(status.stopReason)
  status.completed = true;
  status.stopReason = 'reached_Vmax';
  status.lastVoltage = prevV;
end
if nValid == 0
  status.stopReason = 'empty_curve';
end
end

function [r, ok, failReason] = stepToConverged(dev, Gprofile, prevState, prevV, targetV, depthLeft, JscMax, JscMaxMargin)
opts = struct('psiInit',prevState.psi,'nInit',prevState.n,'pInit',prevState.p);
r = gummelSolveML(dev, targetV, Gprofile, opts);
physicallyValid = ~(abs(targetV) < 1e-12 && abs(r.J) > JscMax*JscMaxMargin);
if r.isFinite && r.converged && physicallyValid
  ok = true; failReason = '';
  return;
end
if ~r.isFinite
  failReason = 'nonfinite';
elseif ~physicallyValid
  failReason = 'current_nonconservation';
else
  failReason = 'nonconvergence';
end
if depthLeft <= 0
  ok = false;
  return;
end
Vmid = 0.5*(prevV + targetV);
[rMid, okMid, failReason] = stepToConverged(dev, Gprofile, prevState, prevV, Vmid, depthLeft-1, JscMax, JscMaxMargin);
if ~okMid
  ok = false;
  return;
end
[r, ok, failReason] = stepToConverged(dev, Gprofile, rMid, Vmid, targetV, depthLeft-1, JscMax, JscMaxMargin);
end
