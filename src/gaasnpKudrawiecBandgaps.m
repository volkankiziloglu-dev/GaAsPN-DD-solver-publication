function out = gaasnpKudrawiecBandgaps(x, pFraction, asFraction, params)
%GAASNPKUDRAWIECBANDGAPS Kudrawiec Gamma/X/L comparator for GaNAsP.
%
% Implements the endpoint-BAC interpolation of Kudrawiec, J. Appl. Phys.
% 101, 116101 (2007), Eq. (3), and the outside-Gamma extension described in
% J. Appl. Phys. 105, 063529 (2009), Eqs. (7)-(8).
%
% Inputs are total group-V mole fractions:
%   x          substitutional nitrogen fraction
%   pFraction  phosphorus fraction
%   asFraction arsenic fraction
%
% The 2009 paper states that X/L coupling is about 3-4 times weaker than
% Gamma.  The published calculation uses a factor of three; this is the
% default here and is exposed in params.outsideGammaReduction.
%
% This function is a validity comparator, not a replacement for the
% project's legacy calibrated BAC expression.

def = struct( ...
  'EgGaAs_Gamma',1.424, 'EgGaP_Gamma',2.78, ...
  'EgGaAs_X',1.981,     'EgGaP_X',2.35, ...
  'EgGaAs_L',1.815,     'EgGaP_L',2.72, ...
  'bowGamma',0.19, 'bowX',0.24, 'bowL',0.16, ...
  'EN_GaAsN',1.65, 'EN_GaPN',2.18, ...
  'C_GaAsN_Gamma',2.70, 'C_GaPN_Gamma',3.05, ...
  'outsideGammaReduction',3.0);
if nargin < 4 || isempty(params), params = struct(); end
fn = fieldnames(def);
for k = 1:numel(fn)
  if ~isfield(params,fn{k}) || isempty(params.(fn{k}))
    params.(fn{k}) = def.(fn{k});
  end
end
p = params;

if any(x(:) < 0 | pFraction(:) < 0 | asFraction(:) < 0)
  error('gaasnpKudrawiecBandgaps:NegativeComposition', ...
    'N, P, and As mole fractions must be non-negative.');
end
sumFraction = x + pFraction + asFraction;
if any(abs(sumFraction(:)-1) > 2e-6)
  error('gaasnpKudrawiecBandgaps:CompositionSum', ...
    'N+P+As must equal one (maximum residual %.3g).', ...
    max(abs(sumFraction(:)-1)));
end
if any(x(:) >= 1)
  error('gaasnpKudrawiecBandgaps:NRange','N fraction must be below one.');
end

% Kudrawiec's z' correction: redistribute N proportionally over the
% non-nitrogen group-V sublattice.  Thus z'=P/(1-N).
zPrime = pFraction ./ (1-x);

CGaAsOutside = p.C_GaAsN_Gamma / p.outsideGammaReduction;
CGaPOutside  = p.C_GaPN_Gamma  / p.outsideGammaReduction;

EgGamma = mixedValley(p.EgGaAs_Gamma,p.EgGaP_Gamma, ...
  p.EN_GaAsN,p.EN_GaPN,p.C_GaAsN_Gamma,p.C_GaPN_Gamma, ...
  p.bowGamma,x,zPrime);
EgX = mixedValley(p.EgGaAs_X,p.EgGaP_X, ...
  p.EN_GaAsN,p.EN_GaPN,CGaAsOutside,CGaPOutside, ...
  p.bowX,x,zPrime);
EgL = mixedValley(p.EgGaAs_L,p.EgGaP_L, ...
  p.EN_GaAsN,p.EN_GaPN,CGaAsOutside,CGaPOutside, ...
  p.bowL,x,zPrime);

allGaps = [EgGamma(:),EgX(:),EgL(:)];
[fundamentalGap,idx] = min(allGaps,[],2);
names = ["Gamma","X","L"];
valley = names(idx);
directMargin = min(EgX,EgL)-EgGamma;
isDirect = directMargin >= 0;

out = struct('EgGamma',EgGamma,'EgX',EgX,'EgL',EgL, ...
  'fundamentalGap',reshape(fundamentalGap,size(EgGamma)), ...
  'fundamentalValley',reshape(valley,size(EgGamma)), ...
  'directGapMargin_eV',directMargin,'isDirect',isDirect, ...
  'zPrime_P_on_nonN_sublattice',zPrime, ...
  'outsideGammaReduction',p.outsideGammaReduction, ...
  'source',"Kudrawiec 2007 DOI 10.1063/1.2736618; " + ...
           "Kudrawiec 2009 DOI 10.1063/1.3087781");
end

function Eg = mixedValley(EAs,EP,ENAs,ENP,CAs,CP,bow,x,zPrime)
EgAsN = lowerBAC(EAs,ENAs,CAs,x);
EgPN  = lowerBAC(EP,ENP,CP,x);
Eg = (1-zPrime).*EgAsN + zPrime.*EgPN - zPrime.*(1-zPrime).*bow;
end

function Eminus = lowerBAC(Ehost,EN,C,x)
Eminus = 0.5*(EN+Ehost-sqrt((EN-Ehost).^2+4*C.^2.*x));
end
