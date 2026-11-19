function [Gnode, transFlux] = singleWavelengthGenML(dev, lam, irr, frontR, cvIntegrate)
% cvIntegrate (optional, default false; 2026-11-07 referee check R1-5): if true, Gnode(i) is the exact
%   Beer-Lambert absorption integrated over the control volume of node i, divided by its width
%   [x_i - h_{i-1}/2, x_i + h_i/2], using the same piecewise-constant attenuation (alpha of the left node
%   of each segment) as the flux recursion; the sum over control volumes then equals the absorbed flux.
%   false reproduces the nodal sampling Gnode(i) = alpha(i)*phi(x_i) of the production runs.
if nargin < 5 || isempty(cvIntegrate), cvIntegrate = false; end
% common/singleWavelengthGen.m'nin cok-katmanli hali: TEK bir isAbs/ETL
% ayrimi yerine, HER dugum kendi soğurma modelini tasir (dev.useOpticalTable
% -> o dugumun KENDI katmanina ait tablosu, dev.opticalTableNode{i};
% yoksa dev.alphaMaxNode/edgeWidthNode ile lojistik-kenar; ikisi de
% yoksa/alphaMax=0 ise seffaf), Beer-Lambert TUM cihaz boyunca kesintisiz
% birikir. 2-bolgeli ozel durumda (ETL seffaf, Abs tabloyla soguruyor)
% common/singleWavelengthGen.m ile AYNI sonucu verir.
%
% frontR: (OPSIYONEL) on-yuzey yansima orani. Ya SKALER (dalga boyundan
%   bagimsiz, eski davranis) ya da bir FONKSIYON TUTAMACI @(lam)->R
%   (dalga-boyu-bagimli, Gorev #10 -- bkz. RearHJT/arcReflectanceGaP.m)
%   olabilir. Bos/verilmemisse: eger 1. dugumun tablosu (veya, geriye-
%   uyumluluk icin, eski tek dev.alphaTable) varsa dogrudan o tablonun
%   ciplak-yuzey Fresnel yansimasi kullanilir (ESKI davranis) -- ONEMLI
%   NOT: bu, cok-katmanli bir yapida (orn. GaP pencere/GaAsNP absorber)
%   fiziksel olarak YANLIS olabilir, cunku isik ONCE pencereye carpar,
%   dogrudan absorbere degil. Dogru sonuc icin cagiran kod frontR'yi
%   ACIKCA (skaler veya fonksiyon) vermelidir -- bkz. RearHJT/run_4layer_with_ARC.m.
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
  frontRVal = ((nFront-1)/(nFront+1))^2; % ciplak-tablo-katmani Fresnel (eski davranis, bkz. yukaridaki not)
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
    a = 0; % bu dugum seffaf
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
