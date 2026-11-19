function p = ganp21Params()
%GANP21PARAMS Nominal binary GaPN layer lattice matched to Si at 300 K.
%
% Optical/material anchor:
% Ben Saddik et al., Mater. Sci. Semicond. Process. (2024),
% DOI 10.1016/j.mssp.2024.109011.
%
% Eg=1.95 eV is an experimental room-temperature anchor.  Electrical
% transport, affinity, doping activation, and lifetime remain uncertain and
% must still be swept/calibrated by each device case.
xN=0.021;
aGaN=4.520;   % angstrom, zincblende GaN
aGaP=5.4505;  % angstrom
aSi=5.4310;   % angstrom
aGaPN=xN*aGaN+(1-xN)*aGaP;

p=struct('N_fraction',xN,'P_fraction',1-xN,'As_fraction',0, ...
  'N_percent',100*xN,'P_percent',100*(1-xN),'As_percent',0, ...
  'Eg',1.95,'a_A',aGaPN,'mismatchToSi_percent',100*(aGaPN-aSi)/aSi, ...
  'source',"Ben Saddik et al. 2024 DOI 10.1016/j.mssp.2024.109011");
end
