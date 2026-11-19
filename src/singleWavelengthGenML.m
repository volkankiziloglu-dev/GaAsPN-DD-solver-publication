function [Gnode, transFlux] = singleWavelengthGenML(dev, lam, irr, frontR, cvIntegrate)
% cvIntegrate (optional, default false; 2026-11-07): if true, Gnode(i) is the exact
%   Beer-Lambert absorption integrated over the control volume of node i, divided by its width
%   [x_i - h_{i-1}/2, x_i + h_i/2], using the same piecewise-constant attenuation (alpha of the left node
%   of each segment) as the flux recursion; the sum over control volumes then equals the absorbed flux.
%   false reproduces the nodal sampling Gnode(i) = alpha(i)*phi(x_i) of the production runs.
if nargin < 5 || isempty(cvIntegrate), cvIntegrate = false; end
% Multilayer version of common/singleWavelengthGen.m: instead of a SINGLE isAbs/ETL
% distinction, EACH node carries its own absorption model (dev.useOpticalTable
% -> the table belonging to that node's OWN layer, dev.opticalTableNode{i};
% otherwise logistic edge via dev.alphaMaxNode/edgeWidthNode; if neither
% exists/alphaMax=0, transparent), Beer-Lambert accumulates continuously through the ENTIRE device.
% In the 2-region special case (ETL transparent, Abs absorbing via table) it gives
% the SAME result as common/singleWavelengthGen.m.
%
% frontR: (OPTIONAL) front-surface reflectance. Can be either a SCALAR (wavelength-
%   independent, old behavior) or a FUNCTION HANDLE @(lam)->R
%   (wavelength-dependent, Task #10 -- see RearHJT/arcReflectanceGaP.m).
%   If empty/not given: if the table of node 1 (or, for backward
%   compatibility, the old single dev.alphaTable) exists, the bare-surface Fresnel
%   reflectance of that table is used directly (OLD behavior) -- IMPORTANT
%   NOTE: in a multilayer structure (e.g. GaP window/GaAsNP absorber) this
%   may be physically WRONG, because light hits the window FIRST,
%   not the absorber directly. For a correct result the calling code must give frontR
%   EXPLICITLY (scalar or function) -- see RearHJT/run_4layer_with_ARC.m.
c = constants();
N = dev.N; x = dev.x; Eg = dev.Eg;
Gnode = zeros(N,1);
if irr <= 0
  transFlux = 0;
  return;
end
Ephoton = c.HC_EVNM / lam;

% Per-NODE optical table lookup (dev.opticalTableNode{i}, from each
% node's OWN layer -- buildDeviceML.m), NOT the old single device-level
% dev.alphaTable. The old mechanism applied ONE shared table to every
% useOpticalTable=true node regardless of layer, silently wrong for any
% device mixing two different real materials (confirmed this session: a
% GaP window intended to use its own real optical data ended up using
% the absorber's GaAsNP table instead, since both had
% useOpticalTable=true and only dev.alphaTable existed). Backward
% compatible: a caller that still sets a single dev.alphaTable and gives
% every table-layer that SAME struct (the historical single-table usage
% pattern, still the common case) gets bit-identical results, since
% dev.opticalTableNode{i} then equals dev.alphaTable at every such node,
% and lookupTableAN below reproduces the original asymmetric boundary
% handling (alpha=0 below Emin, clamped to the table's own last value
% above Emax) exactly.
hasNodeTables = isfield(dev,'opticalTableNode');
hasDevTable = isfield(dev,'alphaTable') && ~isempty(dev.alphaTable);

frontTbl = [];
if hasNodeTables && dev.useOpticalTable(1) && ~isempty(dev.opticalTableNode{1})
  frontTbl = dev.opticalTableNode{1};
elseif hasDevTable
  frontTbl = dev.alphaTable; % legacy fallback for devices built before opticalTableNode existed
end

if nargin >= 4 && ~isempty(frontR) && isa(frontR,'function_handle')
  frontRVal = frontR(lam);
elseif nargin >= 4 && ~isempty(frontR)
  frontRVal = frontR;
elseif ~isempty(frontTbl)
  [~, nFront] = lookupTableAN(frontTbl, Ephoton);
  frontRVal = ((nFront-1)/(nFront+1))^2; % bare-table-layer Fresnel (old behavior, see note above)
else
  frontRVal = 0;
end

phi0 = photonFlux(lam, irr) * (1 - frontRVal);
absDepth = 0;
aNode = zeros(N,1); phiNode = zeros(N,1);
for i = 1:N
  if dev.useOpticalTable(i) && hasNodeTables && ~isempty(dev.opticalTableNode{i})
    a = lookupTableAN(dev.opticalTableNode{i}, Ephoton);
  elseif dev.useOpticalTable(i) && hasDevTable
    a = lookupTableAN(dev.alphaTable, Ephoton); % legacy fallback
  elseif dev.alphaMaxNode(i) > 0
    a = alphaAbsorber(Ephoton, Eg(i), dev.alphaMaxNode(i), dev.edgeWidthNode(i));
  else
    a = 0; % this node is transparent
  end
  fluxAtNode = phi0 * exp(-absDepth);
  Gnode(i) = a * fluxAtNode;
  aNode(i) = a; phiNode(i) = fluxAtNode;
  if i < N
    absDepth = absDepth + a * (x(i+1)-x(i));
  end
end
transFlux = phi0 * exp(-absDepth);
if cvIntegrate
  h = diff(x(:));                                   % segment widths, N-1
  aS = aNode(1:N-1);                                % attenuation of segment k = alpha of its left node
  absHalfR = phiNode(1:N-1) .* (1 - exp(-aS.*h/2));                     % absorbed in [x_k, x_k + h_k/2]
  absHalfL = phiNode(1:N-1) .* exp(-aS.*h/2) .* (1 - exp(-aS.*h/2));    % absorbed in [x_k + h_k/2, x_{k+1}]
  Gnode = zeros(N,1); w = zeros(N,1);
  Gnode(1:N-1) = Gnode(1:N-1) + absHalfR; w(1:N-1) = w(1:N-1) + h/2;
  Gnode(2:N)   = Gnode(2:N)   + absHalfL; w(2:N)   = w(2:N)   + h/2;
  Gnode = Gnode ./ w;
end
end

function [a, n] = lookupTableAN(tbl, Ephoton)
% Reproduces the ORIGINAL (pre-per-node-fix) asymmetric boundary
% handling exactly: alpha=0 below Emin (no measured absorption there,
% treat as transparent), clamped to the table's own last tabulated value
% above Emax (still assumed strongly absorbing past the measured range).
if Ephoton < tbl.Emin
  a = 0; n = tbl.n(1);
elseif Ephoton > tbl.Emax
  a = tbl.alpha(end); n = tbl.n(end);
else
  a = interp1(tbl.E, tbl.alpha, Ephoton, 'linear');
  n = interp1(tbl.E, tbl.n, Ephoton, 'linear');
end
end
