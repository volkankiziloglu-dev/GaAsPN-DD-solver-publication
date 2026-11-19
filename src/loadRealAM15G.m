function pts = loadRealAM15G()
% Loads the real ASTM G173 AM1.5G Global Tilt spectrum from data/AM15.xls.
% Returns Nx2 [lambda_nm, irradiance_Wm2nm].
here = fileparts(mfilename('fullpath'));
persistent cached
if isempty(cached)
  num = xlsread(resolveDataFile('AM15.xls'));
  cached = num(:, [1 2]);
end
pts = cached;
end
