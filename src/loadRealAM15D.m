function pts = loadRealAM15D()
% Loads the real ASTM G173-03 AM1.5 Direct+circumsolar spectrum from
% data/AM15D.csv. Returns Nx2 [lambda_nm, irradiance_Wm2nm].
%
% Source: pvlib-python's bundled ASTMG173.csv (same SMARTS v2.9.2 origin
% as loadRealAM15G.m's data/AM15.xls -- cross-checked this session, the
% "global" column matches our existing AM1.5G data to within numerical
% noise, 0.0000% max relative diff, confirming a shared source). Native
% integrated power 900.14 W/m2 (matches the ASTM standard's 900.1 W/m2),
% NOT rescaled to 1000 W/m2 here -- callers wanting Klitzke 2025's
% calibrated "AM1.5d, 1-sun, 1000 W/m2" measurement condition must scale
% explicitly (spectrumPower-based rescale), same pattern already used for
% the CIE indoor spectra elsewhere in this project.
here = fileparts(mfilename('fullpath'));
persistent cached
if isempty(cached)
  % readmatrix (not readtable): data/AM15D.csv here is headerless
  % (2 columns, no header row) -- readtable would silently misparse the
  % first data row as column names, corrupting the first wavelength
  % point. readmatrix has no such assumption.
  cached = readmatrix(resolveDataFile('AM15D.csv'));
end
pts = cached;
end
