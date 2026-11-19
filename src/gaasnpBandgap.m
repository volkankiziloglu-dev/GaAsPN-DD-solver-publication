function out = gaasnpBandgap(x, params)
% Band-anti-crossing (BAC) bandgap model for GaAs(1-x-y)NxPy lattice-matched to Si.
% Constants and weighting scheme taken verbatim from the user's own GaASNP.m
% (lines 169-177): y,(1-y) used directly as GaPN/GaAsN weights (not renormalized).
%
% x: nitrogen mole fraction (e.g. 0.025-0.05)

def = struct( ...
  'aGaAs', 0.5632, 'aGaP', 0.54505, 'aGaN', 0.4520, 'aSi', 0.5431, ... % nm
  'Eg_GaAs', 1.424, 'Eg_GaP', 2.77, ...   % eV
  'Cb', 0.19, ...                        % eV, GaAsP bowing parameter
  'EN_GaAsN', 1.65, 'EN_GaPN', 2.25, ...  % eV, BAC nitrogen level per host
  'Cmn_GaAsN', 2.70, 'Cmn_GaPN', 3.16 ... % eV, BAC coupling matrix element per host
);
if nargin < 2 || isempty(params), params = struct(); end
fn = fieldnames(def);
for i = 1:numel(fn)
  if ~isfield(params, fn{i}) || isempty(params.(fn{i}))
    params.(fn{i}) = def.(fn{i});
  end
end
p = params;

y = (x.*(p.aGaAs - p.aGaN) + p.aSi - p.aGaAs) ./ (p.aGaP - p.aGaAs);
z = 1 - x - y;

Eg_GaAsP = p.Eg_GaAs.*z + p.Eg_GaP.*(1-z) - z.*(1-z).*p.Cb;

En = y.*p.EN_GaPN + (1-y).*p.EN_GaAsN;
Cmn = y.*p.Cmn_GaPN + (1-y).*p.Cmn_GaAsN;

Eg = 0.5 * (Eg_GaAsP + En - sqrt((Eg_GaAsP - En).^2 + 4*Cmn.^2.*x));

% Independent Kudrawiec (2007/2009) comparator.  This does not replace the
% legacy BAC value above, so existing calibrated/regression cases retain
% their numerical definition.  It supplies composition metadata and a
% Gamma/X/L validity check for the optical direct-gap assumption.
kud = gaasnpKudrawiecBandgaps(x, y, z);

out = struct('x',x,'y',y,'z',z, ...
  'N_fraction',x,'P_fraction',y,'As_fraction',z, ...
  'N_percent',100*x,'P_percent',100*y,'As_percent',100*z, ...
  'composition',compose("GaN%.5fP%.5fAs%.5f",x,y,z), ...
  'Eg_GaAsP',Eg_GaAsP,'EN',En,'Cmn',Cmn,'Eg',Eg, ...
  'kudrawiec',kud,'EgGamma_Kudrawiec',kud.EgGamma, ...
  'EgX_Kudrawiec',kud.EgX,'EgL_Kudrawiec',kud.EgL, ...
  'fundamentalValley_Kudrawiec',kud.fundamentalValley, ...
  'directGapMargin_eV',kud.directGapMargin_eV);
end
