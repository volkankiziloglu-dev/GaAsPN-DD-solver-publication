function [G, transmittedPts] = deviceGenAndTransmitML(dev, spectrumPts, frontR, cvIntegrate)
% cvIntegrate (optional, default false): control-volume-integrated generation, see singleWavelengthGenML.m.
if nargin < 4 || isempty(cvIntegrate), cvIntegrate = false; end
% common/deviceGenAndTransmit.m'nin singleWavelengthGenML.m kullanan hali.
c = constants();
N = dev.N;
G = zeros(N,1);
nLam = size(spectrumPts,1);
transIrr = zeros(nLam,1);

for k = 1:nLam
  lam = spectrumPts(k,1); irr = spectrumPts(k,2);
  [~, tFlux] = singleWavelengthGenML(dev, lam, irr, frontR);
  Ephoton_J = (c.HC_EVNM/lam) * c.Q;
  transIrr(k) = tFlux * Ephoton_J * 1e4;
end
transmittedPts = [spectrumPts(:,1), transIrr];

for k = 1:nLam-1
  l1 = spectrumPts(k,1); v1 = spectrumPts(k,2);
  l2 = spectrumPts(k+1,1); v2 = spectrumPts(k+1,2);
  dl = l2 - l1;
  g1 = singleWavelengthGenML(dev, l1, v1, frontR, cvIntegrate);
  g2 = singleWavelengthGenML(dev, l2, v2, frontR, cvIntegrate);
  G = G + 0.5*(g1+g2)*dl;
end
end
