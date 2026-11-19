function [psi, converged, iterOut, residualNorm] = newtonPoissonML(dev, psiInit, un, vp, psiLeftBC, psiRightBC, equilibriumFlag, verbose)
% Copy of RearHJT/newtonPoissonML.m -- epsL/epsR are now multiplied by c.EPS0
% (see eps0fix_v2/newtonPoisson_v2.m and MATHEMATICAL_MODEL.md
% Section K item 12), PLUS [converged, iterOut, residualNorm] outputs
% (external review item 1.2 -- previously ONLY psi was returned, so the calling code
% could not see whether the inner Poisson convergence actually succeeded).
%
% FIX (external review item 2.2): the convergence decision PREVIOUSLY relied ONLY
% on max(abs(delta))<deltaTol -- a small Newton step does NOT ALWAYS
% mean a small residual (ill-conditioned Jacobian, step
% limiting, stagnation). Now a SCALED residual
% (F normalized by the local charge-density scale) is also
% checked.
if nargin < 8, verbose = false; end
c = constants();
N = dev.N; h = dev.h; eps_ = dev.eps; Nd = dev.Nd; Na = dev.Na;
psi = psiInit;
psi(1) = psiLeftBC; psi(N) = psiRightBC;

idx = 2:N-1;
epsL = c.EPS0 * 0.5*(eps_(idx) + eps_(idx-1));
epsR = c.EPS0 * 0.5*(eps_(idx) + eps_(idx+1));
hL = h(idx-1); hR = h(idx);
hAvg = 0.5*(hL+hR);
coeffL = epsL./(hL.*hAvg);
coeffR = epsR./(hR.*hAvg);
Ndm = Nd(idx); Nam = Na(idx);
offN = dev.psiOffsetN(idx); offP = dev.psiOffsetP(idx);
Nc = dev.Nc(idx); Nv = dev.Nv(idx);
chiRel = dev.chi(idx)-dev.chiRef; Eg = dev.Eg(idx);

a = zeros(N,1); b = zeros(N,1); cc = zeros(N,1); d = zeros(N,1);
a(1)=0; b(1)=1; cc(1)=0; d(1)=0;
a(N)=0; b(N)=1; cc(N)=0; d(N)=0;
a(idx) = coeffL; cc(idx) = coeffR;

converged = false;
deltaTol = 1e-10;
residualTol = 1e-6; % scaled (dimensionless) residual threshold
maxIter = 60;
for iter = 1:maxIter
  if equilibriumFlag
    n_i = Nc .* exp((psi(idx)+chiRel)/c.VT);
    p_i = Nv .* exp(-(psi(idx)+chiRel+Eg)/c.VT);
  else
    psiEffN = psi(idx) + offN;
    psiEffP = psi(idx) + offP;
    n_i = un(idx) .* exp(psiEffN/c.VT);
    p_i = vp(idx) .* exp(-psiEffP/c.VT);
  end
  rho = c.Q * (p_i - n_i + Ndm - Nam);
  F = (epsR.*(psi(idx+1)-psi(idx))./hR - epsL.*(psi(idx)-psi(idx-1))./hL)./hAvg + rho;

  dRhodPsi = c.Q * (-p_i - n_i) / c.VT;

  b(idx) = -(coeffL+coeffR) + dRhodPsi;
  d(idx) = -F;

  delta = thomasSolve(a, b, cc, d);
  % DAMPED NEWTON (new): the inner loop had no relaxation/limiting at all -- after the eps0
  % fix (because of the wide dynamic range of the coefficients),
  % the first Newton step could sometimes make an absurd jump of tens of VOLTS
  % (e.g. 36.75V, see debug_ML_poisson_verbose.m), which led to exp() overflow
  % and NaN. Limiting the step size with maxStep (preserving the direction,
  % scaling down proportionally) is a standard "damped Newton" technique
  % -- it has no effect once converging (delta is already small).
  maxStep = 0.5; % V
  maxAbsDelta = max(abs(delta));
  if maxAbsDelta > maxStep
    delta = delta * (maxStep/maxAbsDelta);
  end
  psi = psi + delta;
  if verbose
    [maxD, iMax] = max(abs(delta));
    fprintf('  newtonPoissonML_v2 iter=%2d  max|delta|=%.4e  at node %d (psi=%.4f)\n', iter, maxD, iMax+1, psi(iMax+1));
  end
  chargeScale = c.Q * max([max(Ndm), max(Nam), max(n_i), max(p_i), 1e10]);
  scaledResidual = max(abs(F)) / chargeScale;
  if max(abs(delta)) < deltaTol && scaledResidual < residualTol
    converged = true;
    break;
  end
end
iterOut = iter;

% re-evaluate F at the final psi and report the true residual (the F BEFORE
% applying delta would show how much correction the last step needed
% -- the residual of the applied psi itself is more accurate here).
if equilibriumFlag
  n_i = Nc .* exp((psi(idx)+chiRel)/c.VT);
  p_i = Nv .* exp(-(psi(idx)+chiRel+Eg)/c.VT);
else
  psiEffN = psi(idx) + offN;
  psiEffP = psi(idx) + offP;
  n_i = un(idx) .* exp(psiEffN/c.VT);
  p_i = vp(idx) .* exp(-psiEffP/c.VT);
end
rhoFinal = c.Q * (p_i - n_i + Ndm - Nam);
Ffinal = (epsR.*(psi(idx+1)-psi(idx))./hR - epsL.*(psi(idx)-psi(idx-1))./hL)./hAvg + rhoFinal;
residualNorm = max(abs(Ffinal));
end
