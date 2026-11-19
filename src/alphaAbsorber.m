function a = alphaAbsorber(E_eV, Eg, alphaMax, width)
% Smoothed absorption edge (logistic step) in cm^-1
a = alphaMax ./ (1 + exp(-(E_eV - Eg)/width));
end
