function [Brad, detail] = computeBradVRS(opt, Eg_eV, Nc_cm3, Nv_cm3, T_K)
%COMPUTEBRADVRS Radiative coefficient from van Roosbroeck--Shockley.
% opt.E [eV], opt.alpha [cm^-1], opt.n [-]; output Brad [cm^3 s^-1].

arguments
  opt (1,1) struct
  Eg_eV (1,1) double {mustBePositive}
  Nc_cm3 (1,1) double {mustBePositive}
  Nv_cm3 (1,1) double {mustBePositive}
  T_K (1,1) double {mustBePositive} = 300
end

kB_eV = 8.617333262145e-5;
h_eVs = 4.135667696e-15;
c_cmps = 2.99792458e10;
VT_eV = kB_eV*T_K;

E = opt.E(:);
alpha = max(opt.alpha(:),0);
nref = opt.n(:);
valid = isfinite(E) & isfinite(alpha) & isfinite(nref) & ...
        E > 0 & alpha >= 0 & nref > 0;
E = E(valid); alpha = alpha(valid); nref = nref(valid);
[E,order] = sort(E);
alpha = alpha(order); nref = nref(order);

% Equilibrium spontaneous-emission rate density integrated over photon
% energy. expm1 avoids cancellation at small E/kT.
rsp_per_eV = alpha .* (2*pi*nref.^2.*E.^2) ./ ...
  (h_eVs^3*c_cmps^2.*expm1(E/VT_eV));
Rrad0_cm3ps = trapz(E,rsp_per_eV);
ni_cm3 = sqrt(Nc_cm3*Nv_cm3)*exp(-Eg_eV/(2*VT_eV));
Brad = Rrad0_cm3ps/ni_cm3^2;

detail = struct('T_K',T_K,'Eg_eV',Eg_eV,'Nc_cm3',Nc_cm3, ...
  'Nv_cm3',Nv_cm3,'ni_cm3',ni_cm3,'Rrad0_cm3ps',Rrad0_cm3ps, ...
  'Emin_eV',min(E),'Emax_eV',max(E),'nSamples',numel(E));
end
