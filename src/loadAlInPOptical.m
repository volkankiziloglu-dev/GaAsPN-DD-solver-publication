function opt = loadAlInPOptical()
%LOADALINPOPTICAL Literature-based optical table for GaAs-matched AlInP.
%   n(lambda): digitized from the n+-AlInP curve in Fig. 8(b) of
%   Ochoa-Martinez et al., Sol. Energy Mater. Sol. Cells 174 (2018)
%   388-396, doi:10.1016/j.solmat.2017.09.028. That sample was
%   Al0.537In0.463P with n=5.6e18 cm^-3, close to the 8e18 cm^-3 Klitzke
%   window used here.
%
%   alpha(lambda): digitized bulk envelope from Fig. 3 of Cheong et al.,
%   J. Phys. D 48 (2015) 405101, doi:10.1088/0022-3727/48/40/405101.
%   The authors report that alpha was insensitive to dopant type over the
%   investigated range. Values below 380 nm follow the high-absorption
%   envelope and values beyond 560 nm are set to zero because they are
%   below the measured indirect edge and negligible for a 25-nm window.
%
%   This is a traceable FIGURE-DIGITIZED dataset, not publisher raw data.
%   It is intended for optical-loss/ARC sensitivity, with uncertainty
%   assessed against the former ideal-transparent window assumption.
here = fileparts(mfilename('fullpath'));
persistent cached
if isempty(cached)
  T = readtable(resolveDataFile('AlInP_Ochoa2018_Cheong2015_optical.csv'));
  E = 1239.842 ./ T.wavelength_nm;
  [E,ord] = sort(E);
  cached = struct('E',E,'alpha',T.alpha_cm1(ord),'n',T.n(ord), ...
    'Emin',min(E),'Emax',max(E), ...
    'provenance','Ochoa-Martinez 2018 Fig. 8(b) n; Cheong 2015 Fig. 3 alpha; figure digitized');
end
opt = cached;
end
