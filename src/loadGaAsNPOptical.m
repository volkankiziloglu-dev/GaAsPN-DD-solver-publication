function opt = loadGaAsNPOptical(N_percent)
% Loads digitized n(E),k(E) data (nk_data: [lambda_um, n, k]) for the given
% nitrogen content (2.5, 3, 3.5, 4, 4.5, or 5 percent) and returns alpha(E)
% [cm^-1] and n(E) sorted by ascending photon energy E [eV].
map = containers.Map({2.5,3,3.5,4,4.5,5}, ...
  {'nk_25N.mat','nk_3N.mat','nk_35N.mat','nk_4N.mat','nk_45N.mat','nk_5N.mat'});
if ~isKey(map, N_percent)
  error('loadGaAsNPOptical: N_percent must be one of 2.5,3,3.5,4,4.5,5 (got %g)', N_percent);
end
here = fileparts(mfilename('fullpath'));
s = load(resolveDataFile(map(N_percent)));
d = s.nk_data; % [lambda_um, n, k]
lam_um = d(:,1); nn = d(:,2); kk = d(:,3);
E = 1.239842 ./ lam_um; % eV

% alpha = 4*pi*k/lambda ; convert lambda um->cm (1um=1e-4cm) for alpha in cm^-1
alpha = 4*pi*kk ./ (lam_um*1e-4);

[E, ord] = sort(E);
alpha = alpha(ord);
nn = nn(ord);

% Digitized data contains duplicate/near-duplicate energy samples; interp1
% requires strictly unique, increasing sample points, so average duplicates.
[Eu, ~, ic] = unique(E);
alphaU = accumarray(ic, alpha, [], @mean);
nU = accumarray(ic, nn, [], @mean);

opt = struct('E', Eu, 'alpha', alphaU, 'n', nU, 'Emin', min(Eu), 'Emax', max(Eu));
end
