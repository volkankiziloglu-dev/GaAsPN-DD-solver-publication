function [F, n, p, JnSeg, JpSeg] = assembleResidualQFPFD12(dev, psi, phin, phip, V, Gprofile)
%ASSEMBLERESIDUALQFPFD12 Fermi-Dirac F_1/2 statistics (fd12.m) with the consistent degenerate Scharfetter-Gummel flux of
% assembleResidualQFPFermiDiracBCSG.m (2026-11-07). Segment degeneracy factor g = (eta2-eta1)/(ln F(eta2)-ln F(eta1))
% (midpoint value 1/(d ln F/d eta) for |eta2-eta1| < 1e-4); the band-parameter part of the drift argument is not divided by g, so the
% flux vanishes exactly in equilibrium, also across heterojunctions. Recombination as in the Blakemore files (n_i,eff supplied by the caller).
c = constants();
N = dev.N;
h = dev.h;

[n, p, etaN, etaP] = carrierDensitiesQFPFD12(dev, psi, phin, phip);

idxP = (2:N-1)';
epsL = c.EPS0*0.5*(dev.eps(idxP)+dev.eps(idxP-1));
epsR = c.EPS0*0.5*(dev.eps(idxP)+dev.eps(idxP+1));
hL = h(idxP-1); hR = h(idxP); hAvgP = 0.5*(hL+hR);
rho = c.Q*(p(idxP)-n(idxP)+dev.Nd(idxP)-dev.Na(idxP));

Fpsi = zeros(N,1);
Fpsi(idxP) = (epsR.*(psi(idxP+1)-psi(idxP))./hR - epsL.*(psi(idxP)-psi(idxP-1))./hL)./hAvgP + rho;
Fpsi(1) = psi(1) - dev.psiEq(1);
Fpsi(N) = psi(N) - (dev.psiEq(N) + V);

idxJ = (1:N-1)';
gN = segG(etaN, idxJ);  gP = segG(etaP, idxJ);
psiEffN = psi + dev.psiOffsetN;
psiEffP = psi + dev.psiOffsetP;
dPsiN = (psiEffN(idxJ+1)-psiEffN(idxJ))/c.VT;
dPsiP = (psiEffP(idxJ+1)-psiEffP(idxJ))/c.VT;
dlnNc = log(dev.Nc(idxJ+1)./dev.Nc(idxJ));  dlnNv = log(dev.Nv(idxJ+1)./dev.Nv(idxJ));
yN = dlnNc + (dPsiN - dlnNc)./gN;
yP = -dlnNv + (dPsiP + dlnNv)./gP;
muNs = 0.5*(dev.muN(idxJ)+dev.muN(idxJ+1));  muPs = 0.5*(dev.muP(idxJ)+dev.muP(idxJ+1));
JnSeg = c.Q*c.VT*muNs.*gN./h(idxJ) .* (n(idxJ+1).*sgBernoulliExpm1(yN) - n(idxJ).*sgBernoulliExpm1(-yN));
JpSeg = c.Q*c.VT*muPs.*gP./h(idxJ) .* (p(idxJ).*sgBernoulliExpm1(yP) - p(idxJ+1).*sgBernoulliExpm1(-yP));

idxC = (2:N-1)';
hLc = h(idxC-1); hRc = h(idxC); hAvgC = 0.5*(hLc+hRc);
npExcess = n(idxC).*p(idxC)-dev.ni(idxC).^2;   % 2026-11-07: Auger added (as assembleResidualQFP.m)
Rbulk = bulkSRHRate(dev, n, p, idxC) + dev.Brad(idxC).*npExcess + ...
  (dev.CnAuger(idxC).*n(idxC)+dev.CpAuger(idxC).*p(idxC)).*npExcess;
assertIfaceModelSupported(dev, mfilename);   % 'four'/'cross' only in assembleResidualQFP (dd_solver v1.1.1)
Rint = interfaceRecombination(dev, n, p, idxC, hAvgC);
Rtot = Rbulk + Rint;

Fphin = zeros(N,1);
Fphip = zeros(N,1);
Fphin(idxC) = (JnSeg(idxC) - JnSeg(idxC-1)) - c.Q*hAvgC.*(Rtot - Gprofile(idxC));
Fphip(idxC) = (JpSeg(idxC) - JpSeg(idxC-1)) + c.Q*hAvgC.*(Rtot - Gprofile(idxC));

Fphin(1) = phin(1);
Fphip(1) = phip(1);

Sback = dev.Sback;
if dev.Na(N) > 0, pRight = dev.Na(N); else, pRight = dev.ni(N)^2/dev.Nd(N); end
nRight = dev.ni(N)^2/pRight;

if dev.Nd(N) > 0
  Fphin(N) = phin(N) - V;
else
  Fphin(N) = JnSeg(N-1) + c.Q*Sback*(n(N)-nRight);
end
if dev.Na(N) > 0
  Fphip(N) = phip(N) - V;
else
  Fphip(N) = JpSeg(N-1) - c.Q*Sback*(p(N)-pRight);
end

F = zeros(3*N,1);
F(1:3:end) = Fpsi;
F(2:3:end) = Fphin;
F(3:3:end) = Fphip;
end

function g = segG(eta, idxJ)
% segment degeneracy factor for F_1/2: log-mean form (exact zero equilibrium flux), midpoint value for small differences (smooth)
e1 = eta(idxJ); e2 = eta(idxJ+1); d = e2 - e1;
[~, dlm] = fd12(0.5*(e1+e2));
g = 1 ./ dlm;
k = abs(d) >= 1e-4;
if any(k)
  [~,~,l1] = fd12(e1(k)); [~,~,l2] = fd12(e2(k));
  g(k) = d(k) ./ (l2 - l1);
end
end
