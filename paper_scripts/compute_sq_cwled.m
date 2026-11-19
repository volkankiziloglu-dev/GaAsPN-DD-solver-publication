function T = compute_sq_cwled()
%COMPUTE_SQ_CWLED Detailed-balance limit under the digitized CW-LED spectrum of Klitzke et al. scaled to the reported 1000 lx power (3.09 W/m2) and 100 lx (0.31 W/m2);
% same assumptions as compute_sq_limits.m (step-function absorber, 300 K, radiative recombination only, front hemisphere only). Eg = 1.857 (GaAsPN) and 1.90 eV (GaInP).
% Output: results/paper_sq_limits/sq_limits_cwled.csv.
assertCanonicalPath(); warning('off','all');
root = fileparts(fileparts(mfilename('fullpath'))); outDir = fullfile(root,'results','paper_sq_limits'); if ~isfolder(outDir), mkdir(outDir); end
q = 1.602176634e-19; h = 6.62607015e-34; c = 2.99792458e8; kB = 1.380649e-23; Tc = 300; kT = kB*Tc;
Tspd = readtable(resolveDataFile('Klitzke2025_CWLED_SPD_digitized.csv'));
base = [Tspd.lambda_nm, Tspd.relative_irradiance];
rows = {};
for Pt = [3.09 0.31]
  s = base; s(:,2) = s(:,2)*(Pt/spectrumPower(base));
  for Eg = [1.857 1.90]
    EgJ = Eg*q;
    Eint = linspace(EgJ, EgJ+60*kT, 20000);
    J0 = q*(2*pi/(h^3*c^2))*trapz(Eint, Eint.^2./(exp(Eint/kT)-1));
    lam = s(:,1); I = s(:,2); Pin = trapz(lam,I);
    Eph = h*c./(lam*1e-9); phi = I./Eph; use = Eph >= EgJ;
    Jsc = q*trapz(lam(use),phi(use)); Voc = (kB*Tc/q)*log(Jsc/J0+1);
    V = linspace(0,Voc,4000); J = Jsc - J0*(exp(V*q/kT)-1); Pmax = max(V.*J);
    rows(end+1,:) = {Eg, Pin, Jsc*0.1*1e3, Voc*1e3, 100*Pmax/(Jsc*Voc), 100*Pmax/Pin}; %#ok<AGROW>
  end
end
T = cell2table(rows,'VariableNames',{'Eg_eV','Pin_Wm2','Jsc_uAcm2','Voc_mV','FF_pct','PCE_pct'});
writetable(T,fullfile(outDir,'sq_limits_cwled.csv')); disp(T);
end
