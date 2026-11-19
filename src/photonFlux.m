function phi = photonFlux(lam_nm, irr_Wm2nm)
% Photon flux [cm^-2 s^-1 nm^-1] from spectral irradiance [W/m2/nm]
c = constants();
Ephoton_J = (c.HC_EVNM ./ lam_nm) * c.Q;
phi = (irr_Wm2nm * 1e-4) ./ Ephoton_J;
end
