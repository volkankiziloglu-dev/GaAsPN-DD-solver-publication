function [RintN, RintP] = interfaceRecombinationNP(dev, n, p, idx, hAvg)
%INTERFACERECOMBINATIONNP Interface sink of the electron and of the hole continuity equation (per node, volumetric).
% For the production (geometric-mean) model and 'onesided' both outputs are the same array, obtained with one call that is
% identical to interfaceRecombination(dev,n,p,idx,hAvg). For 'four' and 'cross' the sinks are placed carrier-specifically
% (see interfaceRecombination.m); every diagnostic that rebuilds the continuity equations must use this function, otherwise it
% reports the interface recombination itself as a residual (dd_solver v1.1.1).
if ifaceCarrierSpecific(dev)
  RintN = interfaceRecombination(dev, n, p, idx, hAvg, 'n');
  RintP = interfaceRecombination(dev, n, p, idx, hAvg, 'p');
else
  RintN = interfaceRecombination(dev, n, p, idx, hAvg);
  RintP = RintN;
end
end
