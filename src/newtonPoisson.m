function [psi, converged, iterOut, residualNorm] = newtonPoisson(dev, psiInit, un, vp, psiLeftBC, psiRightBC, equilibriumFlag, verbose)
% common/newtonPoisson.m'nin kopyasi -- epsL/epsR artik c.EPS0 (bosluk
% gecirgenligi, 8.85e-14 F/cm) ile carpiliyor (bkz. MATHEMATICAL_MODEL.md
% Bolum K madde 12), ARTI [converged, iterOut, residualNorm] cikislari
% (dis inceleme madde 1.2).
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
residualTol = 1e-6; % dis inceleme madde 2.2: olceklendirilmis residual esigi
maxIter = 60;
for iter = 1:maxIter
  if equilibriumFlag
    n_i = Nc .* exp((psi(idx)+chiRel)/c.VT);
    p_i = Nv .* exp(-(psi(idx)+chiRel+Eg)/c.VT);
  else
    n_i = un(idx) .* exp((psi(idx)+offN)/c.VT);
    p_i = vp(idx) .* exp(-(psi(idx)+offP)/c.VT);
  end
  rho = c.Q * (p_i - n_i + Ndm - Nam);
  F = (epsR.*(psi(idx+1)-psi(idx))./hR - epsL.*(psi(idx)-psi(idx-1))./hL)./hAvg + rho;

  dRhodPsi = c.Q * (-p_i - n_i) / c.VT;

  b(idx) = -(coeffL+coeffR) + dRhodPsi;
  d(idx) = -F;

  delta = thomasSolve(a, b, cc, d);
  % DAMPED NEWTON (bkz. newtonPoissonML_v2.m'deki ayni-isimli not).
  maxStep = 0.5; % V
  maxAbsDelta = max(abs(delta));
  if maxAbsDelta > maxStep
    delta = delta * (maxStep/maxAbsDelta);
  end
  psi = psi + delta;
  if verbose
    [maxD, iMax] = max(abs(delta));
    fprintf('  newtonPoisson_v2 iter=%2d  max|delta|=%.4e  at node %d (psi=%.4f)\n', iter, maxD, iMax+1, psi(iMax+1));
  end
  chargeScale = c.Q * max([max(Ndm), max(Nam), max(n_i), max(p_i), 1e10]);
  scaledResidual = max(abs(F)) / chargeScale;
  if max(abs(delta)) < deltaTol && scaledResidual < residualTol
    converged = true;
    break;
  end
end
iterOut = iter;

if equilibriumFlag
  n_i = Nc .* exp((psi(idx)+chiRel)/c.VT);
  p_i = Nv .* exp(-(psi(idx)+chiRel+Eg)/c.VT);
else
  n_i = un(idx) .* exp((psi(idx)+offN)/c.VT);
  p_i = vp(idx) .* exp(-(psi(idx)+offP)/c.VT);
end
rhoFinal = c.Q * (p_i - n_i + Ndm - Nam);
Ffinal = (epsR.*(psi(idx+1)-psi(idx))./hR - epsL.*(psi(idx)-psi(idx-1))./hL)./hAvg + rhoFinal;
residualNorm = max(abs(Ffinal));
end
