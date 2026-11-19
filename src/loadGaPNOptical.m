function opt = loadGaPNOptical()
%LOADGANPOPTICAL Digitized optical data for GaN0.021P0.979 at 300 K.
% Source: Shokhovets et al., Phys. Rev. B 98, 075205 (2018), Fig. 4(b).
% The absorption tail was digitized from the logarithmic inset; n(E) was
% digitized from the main panel. Values are figure-derived, not an author-
% supplied numerical dataset, and should carry digitization uncertainty.
here=fileparts(mfilename('fullpath'));
T=readtable(resolveDataFile('GaPN_Shockhovets2018_digitized.csv'));
opt=struct('E',T.E_eV,'alpha',T.alpha_cm1,'n',T.n, ...
  'Emin',min(T.E_eV),'Emax',max(T.E_eV), ...
  'source','Shokhovets et al. 2018 Fig. 4(b), digitized', ...
  'compositionN',0.021,'temperature_K',300);
end
