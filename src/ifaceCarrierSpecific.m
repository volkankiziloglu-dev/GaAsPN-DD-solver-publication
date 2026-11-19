function tf = ifaceCarrierSpecific(dev)
%IFACECARRIERSPECIFIC True for the opt-in interface models that place the electron and the hole sink of a heterointerface on
% different nodes ('four', 'cross'; see interfaceRecombination.m). False for the production (geometric-mean) model and 'onesided'.
tf = isfield(dev,'ifaceModel') && any(strcmp(dev.ifaceModel, {'four','cross'}));
end
