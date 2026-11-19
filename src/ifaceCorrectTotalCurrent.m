function Jtot = ifaceCorrectTotalCurrent(dev, Jtot, n, p)
%IFACECORRECTTOTALCURRENT Total segment current Jn+Jp with the interface-recombination jump of the carrier-specific models removed.
% With 'four'/'cross' the electron and the hole sink of a heterointerface sit on different nodes, so Jn+Jp is not constant across
% the junction: adding the two continuity rows gives Jtot(i) - Jtot(i-1) = q*(Sn(i) - Sp(i)) at node i (Sn, Sp: areal sinks,
% Jtot(i): segment between nodes i and i+1). Subtracting q*cumsum(Sn - Sp) restores a spatially constant total current, so that
% the conservation checks, the convergence flag and the terminal current are not biased by the interface recombination current.
% Returns Jtot unchanged (bit-identical) for the production model and 'onesided' (dd_solver v1.1.1).
if ~ifaceCarrierSpecific(dev), return; end
c = constants(); N = dev.N; idx = (2:N-1)'; hAvg = 0.5*(dev.h(idx-1) + dev.h(idx));
[RintN, RintP] = interfaceRecombinationNP(dev, n, p, idx, hAvg);
corr = zeros(N-1,1); corr(idx) = -c.Q*cumsum((RintN - RintP).*hAvg);
Jtot = Jtot + reshape(corr, size(Jtot));
end
