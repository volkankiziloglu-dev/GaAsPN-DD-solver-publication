function [wE, wA] = intRecWeights(dev, iE, iA)
%INTRECWEIGHTS Weights of the interface densities, n_int = n(iE)^wE * n(iA)^wA (same for p and ni).
% Default (production): geometric mean, wE = wA = 0.5, exact zero in equilibrium.
% dev.intRecMode = "lowgap" (round-6 review check): densities of the lower-gap (absorber) side node only, also exact in equilibrium.
wE = 0.5; wA = 0.5;
if isfield(dev,'intRecMode') && dev.intRecMode == "lowgap"
  if dev.Eg(iE) <= dev.Eg(iA), wE = 1; wA = 0; else, wE = 0; wA = 1; end
end
end
