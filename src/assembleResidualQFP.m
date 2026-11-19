function [F, n, p, JnSeg, JpSeg] = assembleResidualQFP(dev, psi, phin, phip, V, Gprofile)
%ASSEMBLERESIDUALQFP Fully-coupled residual F(psi,phin,phip)=0.
%
% Unknown packing: x(3*(i-1)+1)=psi(i), x(3*(i-1)+2)=phin(i),
% x(3*(i-1)+3)=phip(i). See docs/MATHEMATICAL_MODEL.md N.16.6/N.17 for the
% derivation. Reuses sgBernoulli.m, bulkSRHRate.m, interfaceRecombination.m
% UNCHANGED (called with the CURRENT n,p, not a frozen "Prev" iterate --
% this residual is exactly nonlinear, no Gummel-style linearization).
%
% Boundary conditions (i=1 front, i=N back):
%   - Front: ALWAYS Dirichlet phin(1)=0, phip(1)=0 -- psi(1) is pinned at
%     dev.psiEq(1) regardless of V, so this exactly reproduces the
%     equilibrium majority/minority densities at the front contact for
%     any V, matching the old nLeftBC/pLeftBC Dirichlet rows.
%   - Back, MAJORITY carrier: Dirichlet phi(N)=V, NOT phi(N)=0. Derived by
%     requiring n(N) (or p(N)) stay pinned at its FIXED equilibrium value
%     while psi(N)=dev.psiEq(N)+V shifts with bias -- an ideal ohmic
%     contact's quasi-Fermi level moves in lockstep with the local
%     vacuum-level shift the contact itself imposes. (A naive phi(N)=0
%     would leave the majority density V-independent in the WRONG way --
%     it would silently break the applied-bias dependence for that
%     carrier. Verified by hand against the old code's fixed nRightBC/
%     pRightBC Dirichlet values plus solveEquilibrium.m's own phin=phip=0
%     equilibrium definition.)
%   - Back, MINORITY carrier (finite Sback): nonlinear Robin residual,
%     algebraically re-derived from solveElectronContinuity.m's/
%     solveHoleContinuity.m's own back-row coefficients (their comments'
%     stated physical law Jn/q=-Sback*(n-n0)) -- confirmed to collapse to
%     exactly this form by expanding their a(N)/b(N)/d(N) triple.
c = constants();
N = dev.N;
h = dev.h;

[n, p] = carrierDensitiesQFP(dev, psi, phin, phip);

% ---- Poisson (interior nodes) ----
idxP = (2:N-1)';
epsL = c.EPS0*0.5*(dev.eps(idxP)+dev.eps(idxP-1));
epsR = c.EPS0*0.5*(dev.eps(idxP)+dev.eps(idxP+1));
hL = h(idxP-1); hR = h(idxP); hAvgP = 0.5*(hL+hR);
rho = c.Q*(p(idxP)-n(idxP)+dev.Nd(idxP)-dev.Na(idxP));

Fpsi = zeros(N,1);
Fpsi(idxP) = (epsR.*(psi(idxP+1)-psi(idxP))./hR - epsL.*(psi(idxP)-psi(idxP-1))./hL)./hAvgP + rho;
Fpsi(1) = psi(1) - dev.psiEq(1);
Fpsi(N) = psi(N) - (dev.psiEq(N) + V);

% ---- Segment (Scharfetter-Gummel) currents, i=1..N-1 ----
idxJ = (1:N-1)';
psiEffN = psi + dev.psiOffsetN;
psiEffP = psi + dev.psiOffsetP;
dPsiN = (psiEffN(idxJ+1)-psiEffN(idxJ))/c.VT;
dPsiP = (psiEffP(idxJ+1)-psiEffP(idxJ))/c.VT;
DnSeg = c.VT*0.5*(dev.muN(idxJ)+dev.muN(idxJ+1));
DpSeg = c.VT*0.5*(dev.muP(idxJ)+dev.muP(idxJ+1));
JnSeg = c.Q*DnSeg./h(idxJ) .* (n(idxJ+1).*sgBernoulli(dPsiN) - n(idxJ).*sgBernoulli(-dPsiN));
JpSeg = c.Q*DpSeg./h(idxJ) .* (p(idxJ).*sgBernoulli(dPsiP) - p(idxJ+1).*sgBernoulli(-dPsiP));

% ---- Continuity (interior nodes) ----
idxC = (2:N-1)';
hLc = h(idxC-1); hRc = h(idxC); hAvgC = 0.5*(hLc+hRc);
npExcess = n(idxC).*p(idxC)-dev.ni(idxC).^2;
Rbulk = bulkSRHRate(dev, n, p, idxC) + dev.Brad(idxC).*npExcess + ...
  (dev.CnAuger(idxC).*n(idxC)+dev.CpAuger(idxC).*p(idxC)).*npExcess;
if isfield(dev,'ifaceModel') && any(strcmp(dev.ifaceModel,{'cross','four'}))
  % opt-in models with carrier-specific sink placement (dd_solver v1.1.0, see interfaceRecombination.m)
  RtotN = Rbulk + interfaceRecombination(dev, n, p, idxC, hAvgC, 'n');
  RtotP = Rbulk + interfaceRecombination(dev, n, p, idxC, hAvgC, 'p');
else
  Rint = interfaceRecombination(dev, n, p, idxC, hAvgC);
  RtotN = Rbulk + Rint; RtotP = RtotN;
end

% Scaled in CURRENT-DENSITY units (A/cm^2), not the "raw" volumetric-rate
% units (cm^-3 s^-1) div(Jn)/q-(R-G)=0 is usually written in. Multiplying
% the whole equation by Q*hAvg(i) is mathematically exact (same root) but
% keeps this row's natural magnitude comparable to JnSeg/JpSeg themselves
% (~1e-3 to ~10 A/cm^2 for realistic operating points) instead of
% amplifying it by ~1/(Q*hAvg) ~ 1e25-1e26 -- at heavily-doped nodes
% (n~1e18-1e19) that amplification alone was enough to make the combined
% {psi,phin,phip} Jacobian numerically singular (RCOND~1e-35) even though
% every UNKNOWN individually is well-scaled. This is a distinct
% conditioning issue from N.16's original raw-density problem: here the
% UNKNOWNS are fine, but the EQUATIONS themselves spanned ~30 orders of
% magnitude in natural units before this rescaling.
Fphin = zeros(N,1);
Fphip = zeros(N,1);
Fphin(idxC) = (JnSeg(idxC) - JnSeg(idxC-1)) - c.Q*hAvgC.*(RtotN - Gprofile(idxC));
Fphip(idxC) = (JpSeg(idxC) - JpSeg(idxC-1)) + c.Q*hAvgC.*(RtotP - Gprofile(idxC));

% ---- Front boundary (i=1): Dirichlet phi=0 (ideal ohmic contact) ----
% Optional (2026-11-07, referee check): a finite dev.Sfront replaces the
% minority-carrier Dirichlet row by the Robin law -J/q = Sfront*(c - c0),
% the mirror of the back contact; without the field (or Sfront = Inf) the
% residual is unchanged.
Sfront = Inf; if isfield(dev,'Sfront') && ~isempty(dev.Sfront), Sfront = dev.Sfront; end
Fphin(1) = phin(1);
Fphip(1) = phip(1);
if isfinite(Sfront)
  if dev.Na(1) > 0, pLeft = dev.Na(1); else, pLeft = dev.ni(1)^2/dev.Nd(1); end
  nLeft = dev.ni(1)^2/pLeft;
  if dev.Nd(1) <= 0, Fphin(1) = -JnSeg(1) + c.Q*Sfront*(n(1)-nLeft); end
  if dev.Na(1) <= 0, Fphip(1) = -JpSeg(1) - c.Q*Sfront*(p(1)-pLeft); end
end

% ---- Back boundary (i=N) ----
Sback = dev.Sback;
if dev.Na(N) > 0, pRight = dev.Na(N); else, pRight = dev.ni(N)^2/dev.Nd(N); end
nRight = dev.ni(N)^2/pRight;

% Same current-density-unit rescaling as the interior rows above
% (multiply the "Jn/q=-Sback*(n-n0)" law through by Q).
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
