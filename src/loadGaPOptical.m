function opt = loadGaPOptical()
% Loads real GaP n(E),k(E) data (Aspnes & Studna 1983, Phys. Rev. B 27,
% 985-1009, <110> undoped room temperature -- via refractiveindex.info's
% public-domain (CC0) database) and returns alpha(E) [cm^-1] and n(E)
% sorted by ascending photon energy E [eV], same struct format as
% loadGaAsNPOptical.m (E, alpha, n, Emin, Emax) so it is a drop-in
% opticalTable for any GaP window/BSF layer.
%
% Motivation: every #47-#60-family script's GaP window/BSF layers use
% useOpticalTable=false, alphaMax=0 (fully transparent at every photon
% energy) -- confirmed this session to produce model IQE > 1 (physically
% impossible) across most of #58's tested spectrum (2.0-3.55eV, up to
% 1.17), the same defect class already documented for MD421's GaP layers
% above ~3.1eV (MATHEMATICAL_MODEL.md N.16.13/14). This loader supplies
% the missing physical absorption; NOT yet wired into any production
% device script (opt-in, see the caller for exact usage).
here = fileparts(mfilename('fullpath'));
persistent cached
if isempty(cached)
  T = readtable(resolveDataFile('GaP_Aspnes1983_nk.csv'));
  lam_um = T.lambda_um; nn = T.n; kk = T.k;
  E = 1.239842 ./ lam_um; % eV

  alpha = 4*pi*kk ./ (lam_um*1e-4); % cm^-1

  [E, ord] = sort(E);
  alpha = alpha(ord);
  nn = nn(ord);

  [Eu, ~, ic] = unique(E);
  alphaU = accumarray(ic, alpha, [], @mean);
  nU = accumarray(ic, nn, [], @mean);

  cached = struct('E', Eu, 'alpha', alphaU, 'n', nU, 'Emin', min(Eu), 'Emax', max(Eu));
end
opt = cached;
end
