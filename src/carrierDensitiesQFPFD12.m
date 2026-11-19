function [n, p, etaN, etaP] = carrierDensitiesQFPFD12(dev, psi, phin, phip)
%CARRIERDENSITIESQFPFD12 Fermi-Dirac (F_1/2, fd12.m) carrier densities from the quasi-Fermi potentials (same reduced Fermi levels as
% carrierDensitiesQFPFermiDirac.m, which uses the Blakemore approximation): n = Nc F(eta_n), p = Nv F(eta_p).
c = constants();
chiRel = dev.chi - dev.chiRef;
etaN = (psi + chiRel - phin) / c.VT;
etaP = -(psi + chiRel + dev.Eg - phip) / c.VT;
n = dev.Nc .* fd12(etaN);
p = dev.Nv .* fd12(etaP);
end
