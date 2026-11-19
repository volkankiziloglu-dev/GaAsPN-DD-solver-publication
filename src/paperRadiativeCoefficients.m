function B = paperRadiativeCoefficients(T_K)
%PAPERRADIATIVECOEFFICIENTS VRS coefficients for the paper materials.
arguments
  T_K (1,1) double {mustBePositive} = 300
end

absBand = gaasnpBandgap(0.03);
optAbs = loadGaAsNPOptical(3.0);
optGaP = loadGaPOptical();
optGaPN = loadGaPNOptical();

[B.GaAsPN_cm3ps,B.GaAsPN_detail] = computeBradVRS( ...
  optAbs,absBand.Eg,1e19,1e19,T_K);
[B.GaP_cm3ps,B.GaP_detail] = computeBradVRS( ...
  optGaP,2.26,1.8e19,1.9e19,T_K);

[B.GaPN_cm3ps,B.GaPN_detail] = computeBradVRS( ...
  optGaPN,1.95,1.8e19,1.9e19,T_K);
% Diagnostic separation: the measured sub-2.1-eV tail is associated with
% localized/disorder states and is not unambiguously a band-to-band np
% channel. Report cutoff variants rather than hiding this model choice.
opt18=trimOptical(optGaPN,1.8);
opt21=trimOptical(optGaPN,2.1);
B.GaPN_Ege18_cm3ps=computeBradVRS(opt18,1.95,1.8e19,1.9e19,T_K);
B.GaPN_Ege21_cm3ps=computeBradVRS(opt21,1.95,1.8e19,1.9e19,T_K);
B.GaPN_optical_basis = "Shokhovets 2018 Fig. 4(b), digitized";
B.method = "van Roosbroeck--Shockley detailed balance";
B.temperature_K = T_K;
end

function out=trimOptical(opt,Emin)
keep=opt.E>=Emin;
out=struct('E',opt.E(keep),'alpha',opt.alpha(keep),'n',opt.n(keep));
end
