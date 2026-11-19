function T = compute_sq_limits()
%COMPUTE_SQ_LIMITS Detailed-balance (Shockley-Queisser type) limits for the bandgaps and spectra of the paper.
% Step-function absorber (unit absorptance above Eg, none below), 300 K, radiative recombination only, emission into the front hemisphere only (ideal back
% reflector), no reflection, no series/shunt resistance. Spectra: ASTM G173 AM1.5G, and the CIE LED/fluorescent distributions scaled to 1000 lx
% (buildCIEIndoorSpectrum, same as the device calculations). Pin from the full spectrum, as for the device efficiencies.
% Output: results/paper_sq_limits/sq_limits.csv.
assertCanonicalPath(); warning('off','all');
root = fileparts(fileparts(mfilename('fullpath'))); outDir = fullfile(root,'results','paper_sq_limits'); if ~isfolder(outDir), mkdir(outDir); end
q = 1.602176634e-19; h = 6.62607015e-34; c = 2.99792458e8; kB = 1.380649e-23; Tc = 300; kT = kB*Tc;
Egs = [1.857 1.90 1.95];
names = {'AM1.5G','2700K','4000K','6500K','4000F'};
specs = {loadRealAM15G(), buildCIEIndoorSpectrum('2700K',1000), buildCIEIndoorSpectrum('4000K',1000), buildCIEIndoorSpectrum('6500K',1000), buildCIEIndoorSpectrum('4000F',1000)};
rows = {};
for ig = 1:numel(Egs)
  Eg = Egs(ig); EgJ = Eg*q;
  % radiative saturation current density, front hemisphere only: J0 = q*(2*pi/(h^3 c^2)) * int_{Eg}^{inf} E^2/(exp(E/kT)-1) dE   [A/m^2]
  Eint = linspace(EgJ, EgJ+60*kT, 20000);
  J0 = q*(2*pi/(h^3*c^2))*trapz(Eint, Eint.^2./(exp(Eint/kT)-1));
  for is = 1:numel(specs)
    s = specs{is}; lam = s(:,1); I = s(:,2);                  % nm, W m^-2 nm^-1
    Pin = trapz(lam,I);
    Eph = h*c./(lam*1e-9);                                    % J
    phi = I./Eph;                                             % photons m^-2 s^-1 nm^-1
    use = Eph >= EgJ;
    Jsc = q*trapz(lam(use),phi(use));                         % A/m^2
    Voc = (kB*Tc/q)*log(Jsc/J0+1);
    V = linspace(0,Voc,4000); J = Jsc - J0*(exp(V*q/kT)-1);
    [Pmax,im] = max(V.*J);
    FF = Pmax/(Jsc*Voc); PCE = 100*Pmax/Pin;
    rows(end+1,:) = {Eg, names{is}, Pin, Jsc*0.1, Voc*1e3, FF*100, PCE}; %#ok<AGROW>
  end
end
T = cell2table(rows,'VariableNames',{'Eg_eV','Spectrum','Pin_Wm2','Jsc_mAcm2','Voc_mV','FF_pct','PCE_pct'});
writetable(T,fullfile(outDir,'sq_limits.csv')); disp(T);
end
