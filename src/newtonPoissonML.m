function [psi, converged, iterOut, residualNorm] = newtonPoissonML(dev, psiInit, un, vp, psiLeftBC, psiRightBC, equilibriumFlag, verbose)
% RearHJT/newtonPoissonML.m'nin kopyasi -- epsL/epsR artik c.EPS0 ile
% carpiliyor (bkz. eps0fix_v2/newtonPoisson_v2.m ve MATHEMATICAL_MODEL.md
% Bolum K madde 12), ARTI [converged, iterOut, residualNorm] cikislari
% (dis inceleme madde 1.2 -- onceden SADECE psi donuyordu, cagiran kod
% ic Poisson yakinsamasinin gercekte basarili olup olmadigini goremiyordu).
%
% DUZELTME (dis inceleme madde 2.2): yakinsama karari ONCEDEN SADECE
% max(abs(delta))<deltaTol'a dayaniyordu -- kucuk bir Newton adimi HER
% ZAMAN kucuk bir residual anlamina gelmez (kotu kosullu Jacobian, adim
% sinirlamasi, stagnasyon). Simdi ek olarak OLCEKLENDIRILMIS residual
% (F'nin, yerel yuk yogunlugu olcegine gore normlanmis hali) de kontrol
% ediliyor.
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
residualTol = 1e-6; % olceklendirilmis (boyutsuz) residual esigi
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
  % DAMPED NEWTON (yeni): ic dongude hic gevseme/sinirlama yoktu -- eps0
  % duzeltmesinden sonra (katsayilarin genis dinamik araligi yuzunden),
  % ilk Newton adimi bazen onlarca VOLT'luk sacma bir sicrama yapabiliyor
  % (orn. 36.75V, bkz. debug_ML_poisson_verbose.m), bu da exp() tasmasina
  % ve NaN'a yol aciyordu. Adim buyuklugunu maxStep ile sinirlamak (yonu
  % koruyarak, orantili kucultmek), standart bir "damped Newton" teknigi
  % -- yakinsadiginda (delta zaten kucuk) hicbir etkisi yok.
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

% son psi'de F'yi yeniden degerlendirip gercek residual'i rapor et (delta
% uygulanmadan ONCEKI F, son adimin ne kadar duzeltme gerektirdigini
% gosterirdi -- burada uygulanmis psi'nin kendi residual'i daha dogru).
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
