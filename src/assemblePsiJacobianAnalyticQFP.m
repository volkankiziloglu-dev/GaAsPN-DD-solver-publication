function J = assemblePsiJacobianAnalyticQFP(dev, psi, phin, phip)
%ASSEMBLEPSIJACOBIANANALYTICQFP Exact d(F_psi)/d(psi), phin/phip frozen.
% Tridiagonal (interior rows only depend on psi(i-1),psi(i),psi(i+1)).
% Mirrors assembleResidualQFP.m's Fpsi construction exactly -- see that
% file for the residual formula this differentiates. Structurally
% identical to newtonPoissonML.m's existing analytic Poisson Jacobian
% (same physics, same discretization), re-derived here directly from
% carrierDensitiesQFP.m so the two Jacobian paths (FD vs analytic) stay
% independently traceable to the same residual source.
c = constants();
N = dev.N;
h = dev.h;
[n, p] = carrierDensitiesQFP(dev, psi, phin, phip);

chiRel = dev.chi - dev.chiRef;
argN = (psi + chiRel - phin) / c.VT;
argP = -(psi + chiRel + dev.Eg - phip) / c.VT;
EXP_CLAMP = 700;
clampedN = abs(argN) >= EXP_CLAMP;
clampedP = abs(argP) >= EXP_CLAMP;
dn_dpsi = n / c.VT; dn_dpsi(clampedN) = 0;
dp_dpsi = -p / c.VT; dp_dpsi(clampedP) = 0;

idxP = (2:N-1)';
epsL = c.EPS0*0.5*(dev.eps(idxP)+dev.eps(idxP-1));
epsR = c.EPS0*0.5*(dev.eps(idxP)+dev.eps(idxP+1));
hL = h(idxP-1); hR = h(idxP); hAvgP = 0.5*(hL+hR);

lower = zeros(N,1); diagv = zeros(N,1); upper = zeros(N,1);
diagv(1) = 1; % Fpsi(1) = psi(1) - psiEq(1)
diagv(N) = 1; % Fpsi(N) = psi(N) - (psiEq(N)+V)

drho_dpsi = c.Q * (dp_dpsi(idxP) - dn_dpsi(idxP));
lower(idxP) = epsL ./ hL ./ hAvgP;
upper(idxP) = epsR ./ hR ./ hAvgP;
diagv(idxP) = -epsR./hR./hAvgP - epsL./hL./hAvgP + drho_dpsi;

J = sparse(N,N);
rows = (1:N)';
J = J + sparse(rows, rows, diagv, N, N);
J = J + sparse(idxP, idxP-1, lower(idxP), N, N);
J = J + sparse(idxP, idxP+1, upper(idxP), N, N);
end
