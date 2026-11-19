function psi = solveEquilibriumFD12(dev)
%SOLVEEQUILIBRIUMFD12 (2026-11-07) as solveEquilibriumFermiDirac.m with F_1/2 statistics (fd12.m, assembleResidualQFPFD12.m).
%SOLVEEQUILIBRIUMFERMIDIRAC Constant-Fermi equilibrium under Blakemore-
% corrected carrier statistics (carrierDensitiesQFPFermiDirac.m),
% mirroring solveEquilibrium.m's role but with the FD-correct inverse.
%
% Why this is needed (not just cosmetic): solveEquilibrium.m's own
% per-node initial guess inverts n0 -> psi via the BOLTZMANN relation
% psi+chiRel=VT*log(n0/Nc). Reusing that psi as the Dirichlet boundary
% target for a Fermi-Dirac-corrected residual (assembleResidualQFPFermiDirac.m)
% left EFn/EFp badly non-flat at true equilibrium (~5-8mV span measured
% on #58's own degenerate window/rear layers, Nd/Nc~0.556, Na/Nv~0.526,
% see docs/THERMIONIC_INTERFACE_PLAN.md Section 5c / MATHEMATICAL_MODEL.md
% N.18) -- because those same layers are exactly where boundary psi(1)/
% psi(N) are set, and Boltzmann's inversion is measurably wrong there.
%
% Charge neutrality (n(psi)-p(psi)-Dop=0) does NOT depend on which
% carrier-density formula is used for the MAGNITUDE of the majority
% concentration when doping is strongly one-sided (n0~=Nd or p0~=Na to
% leading order, since the minority term is negligible either way) --
% solveEquilibrium.m's own majority/minority quadratic (from n0-p0=Dop,
% n0*p0~=ni^2 for the minority estimate) is kept AS THE INITIAL GUESS
% for n0/p0. The only change is the INVERSION n0 -> psi, corrected to
% the Blakemore form: psi+chiRel = -VT*log(Nc/n0 - gamma) (derived by
% inverting carrierDensitiesQFPFermiDirac.m's own n=Nc/(exp(-eta)+gamma)
% forward formula), instead of Boltzmann's psi+chiRel=VT*log(n0/Nc).
%
% The interior is then relaxed via a psi-ONLY Newton-Poisson solve
% (phin=phip=0 fixed throughout, pure electrostatic equilibrium) using
% assembleResidualQFPFermiDirac.m's own Fpsi block and FD Jacobian
% (assembleJacobianFD.m, unchanged) -- NOT the old Boltzmann-specific
% newtonPoisson.m.
c = constants();
N = dev.N;
gamma = 0.27;
chiRel = dev.chi-dev.chiRef;

psi = zeros(N,1);
for i = 1:N
  Dop = dev.Nd(i)-dev.Na(i);
  majority = 0.5*(abs(Dop)+sqrt(Dop^2+4*dev.ni(i)^2));
  minority = dev.ni(i)^2/majority; %#ok<NASGU>
  if Dop >= 0
    n0 = majority;
    argInv = max(dev.Nc(i)/n0 - gamma, 1e-300);
    psi(i) = c.VT*fd12inv(n0/dev.Nc(i)) - chiRel(i);
  else
    p0 = majority;
    argInv = max(dev.Nv(i)/p0 - gamma, 1e-300);
    % Invert p=Nv/(exp(-etaP)+gamma), etaP=-(psi+chiRel+Eg)/VT:
    % etaP = -log(argInv) -> psi = -chiRel-Eg-VT*etaP = -chiRel-Eg+VT*log(argInv)
    psi(i) = -chiRel(i) - dev.Eg(i) - c.VT*fd12inv(p0/dev.Nv(i));
  end
end

psiTarget1 = psi(1); psiTargetN = psi(N);
phinZero = zeros(N,1); phipZero = zeros(N,1);
Gzero = zeros(N,1);
devTmp = dev; devTmp.psiEq = psi; % Fpsi's boundary rows read dev.psiEq(1)/(N)

resFun = @(psiVec) assembleResidualQFPFD12(devTmp, psiVec, phinZero, phipZero, 0, Gzero);
rowIdx = 1:3:3*N; % psi rows only
maxIter = 60; tol = 1e-13; maxStep = 1.0;
for iter = 1:maxIter
  Ffull = resFun(psi);
  Fpsi = Ffull(rowIdx);
  if ~all(isfinite(Fpsi))
    error('solveEquilibriumFermiDirac:nonfinite', 'Non-finite residual during equilibrium solve.');
  end
  jacFunFull = @(pv) subsetRows(resFun(pv), rowIdx);
  F0forJ = jacFunFull(psi);
  J = assembleJacobianFD(jacFunFull, psi, F0forJ, N, 1, 1);
  delta = J \ (-Fpsi);
  if ~all(isfinite(delta))
    error('solveEquilibriumFermiDirac:nonfiniteDelta', 'Non-finite Newton step during equilibrium solve.');
  end
  delta = max(min(delta, maxStep), -maxStep);
  psi = psi + delta;
  psi(1) = psiTarget1; psi(N) = psiTargetN;
  if max(abs(delta)) < tol && iter > 1
    break;
  end
end
end

function v = subsetRows(Ffull, rowIdx)
v = Ffull(rowIdx);
end
