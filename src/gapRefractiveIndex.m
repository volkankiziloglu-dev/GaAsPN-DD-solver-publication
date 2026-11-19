function n = gapRefractiveIndex(lambda_nm)
% GaP refractive index n(lambda), RELATIVELY ROUGH/APPROXIMATE -- BASED on the
% known general trend of Aspnes & Studna (1983)-type measurements (indirect-gap III-V,
% dispersion increasing toward the 550nm/2.26eV band edge, flattening to ~3.2 at long
% wavelengths), but NOT fitted point by point to a COMPLETE source table
% -- should be treated as a "reasonable starting point" in the SAME spirit as the chi
% uncertainty in RearHJT/gapParams.m. Used only for the ARC
% reflectance calculation (Task #10); GaP still KEEPS the alphaMax=0
% (optically transparent) assumption, this only adds n (for reflectance).
lamTab = [400 450 500 550 600 650 700 750 800 850 900];
nTab   = [3.90 3.75 3.60 3.55 3.50 3.42 3.35 3.29 3.24 3.20 3.16];
n = interp1(lamTab, nTab, lambda_nm, 'linear', 'extrap');
end
