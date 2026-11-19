function R = arcReflectanceGaP(lambda_nm)
% Two-layer ARC (Task #10): MgF2 (outer, air side) / ZnS
% (inner, GaP window side) -- a standard design in the literature for III-V (GaAs) solar cells
% (source: MgF2/ZnS double-layer ARC,
% reported to give R<%10 over 450-1730nm for GaAs single-junction cells
% with the d_MgF2=105nm / d_ZnS=89nm design -- see research notes).
% Since the refractive index of the GaP substrate is of similar magnitude (n~3.2-3.6),
% this design is also a REASONABLE starting point for our GaP-front structure, but it was
% NOT optimized SPECIFICALLY for GaP -- as a further step (the "next step" of
% Task #10) the thicknesses can be re-optimized for GaP's own
% refractive index.
n_air = 1.0;
n_MgF2 = 1.38; d_MgF2 = 105; % nm, outer (air side)
n_ZnS = 2.35;  d_ZnS = 89;   % nm, inner (GaP side)
n_GaP = gapRefractiveIndex(lambda_nm);
R = tmmReflectance(lambda_nm, [n_air, n_MgF2, n_ZnS, n_GaP], [d_MgF2, d_ZnS]);
end
