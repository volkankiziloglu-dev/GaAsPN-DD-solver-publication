function psi = solveEquilibrium(dev)
%SOLVEEQUILIBRIUM Constant-Fermi equilibrium for heterogeneous materials.
c = constants();
N = dev.N;
psi = zeros(N,1);
chiRel = dev.chi-dev.chiRef;

% Gauge: shifted equilibrium Fermi level is zero.
for i = 1:N
  Dop = dev.Nd(i)-dev.Na(i);
  majority = 0.5*(abs(Dop)+sqrt(Dop^2+4*dev.ni(i)^2));
  minority = dev.ni(i)^2/majority;
  if Dop >= 0, n0=majority; else, n0=minority; end
  psi(i) = c.VT*log(n0/dev.Nc(i))-chiRel(i);
end

psi = newtonPoisson(dev,psi,[],[],psi(1),psi(N),true);
end
