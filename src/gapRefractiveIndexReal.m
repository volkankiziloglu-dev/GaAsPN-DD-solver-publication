function n = gapRefractiveIndexReal(lambda_nm)
% Real GaP refractive index n(lambda), from loadGaPOptical.m's own
% Aspnes & Studna 1983 data (src/data/GaP_Aspnes1983_nk.csv,
% 206.6-826.6nm) -- replaces gapRefractiveIndex.m's hand-typed,
% unverified 11-point table (flagged in MATHEMATICAL_MODEL.md L.4 as
% "an unverified, order-of-magnitude-typical estimate, not fit to
% source data") for callers that want the real measured value instead.
% Linear extrapolation beyond 826.6nm (n is slowly, monotonically
% decreasing there in the real data, so a short extrapolation to
% ~900-1000nm is low-risk).
opt = loadGaPOptical();
lam_um = 1.239842 ./ opt.E; % back to wavelength, um
[lam_um, ord] = sort(lam_um);
nSorted = opt.n(ord);
n = interp1(lam_um*1e3, nSorted, lambda_nm, 'linear', 'extrap'); % lam_um*1e3 = nm
end
