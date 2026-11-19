function [n, p] = carrierDensitiesQFP(dev, psi, phin, phip)
%CARRIERDENSITIESQFP Carrier densities from quasi-Fermi potentials.
%
% n(x) = Nc(x)*exp((psi(x)+chiRel(x)-phin(x))/VT)
% p(x) = Nv(x)*exp(-(psi(x)+chiRel(x)+Eg(x)-phip(x))/VT)
%
% chiRel = dev.chi-dev.chiRef. phin=phip=0 reduces this exactly to
% solveEquilibrium.m's own n0/p0 formulas. Cross-checked algebraically
% against gummelSolveML.m's existing un/vp bridge variables
% (psiEffN=psi+dev.psiOffsetN, un=n.*exp(-psiEffN/VT)): un = NcAbs*exp(-phin/VT).
% See docs/MATHEMATICAL_MODEL.md N.16.6/N.17.
c = constants();
EXP_CLAMP = 700; % same defensive clamp style as gummelSolveML.m/newtonPoissonML.m
chiRel = dev.chi - dev.chiRef;
argN = (psi + chiRel - phin) / c.VT;
argP = -(psi + chiRel + dev.Eg - phip) / c.VT;
argN = max(min(argN, EXP_CLAMP), -EXP_CLAMP);
argP = max(min(argP, EXP_CLAMP), -EXP_CLAMP);
n = dev.Nc .* exp(argN);
p = dev.Nv .* exp(argP);
end
