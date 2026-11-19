function dev = buildDeviceML(layers, absorberIdx, junctions, opts)
% General N-layer device builder. Generalized version of buildDevice.m (common/, 2-region
% ETL/Abs) -- as a separate module, WITHOUT TOUCHING the existing
% file. Its output (dev struct) can be used DIRECTLY by solveElectronContinuity.m
% and solveHoleContinuity.m in common/ (same field
% names/formats), because these two files are already written in a way that
% generalizes to any number of regions (per-node arrays only).
%
% layers: struct array, each with:
%   t (thickness, cm), N (number of grid points), Eg, chi, eps, Nd, Na
%   (doping, OPTIONAL for either one, default 0), muN, muP, tau, NcNv,
%   Brad (OPTIONAL, default 0), alphaMax, edgeWidth (OPTIONAL,
%   for the logistic-edge absorption model; if unused, default 0/0.05 --
%   i.e. that layer is optically TRANSPARENT), useOpticalTable (OPTIONAL, if true
%   the real n(E),k(E) table is used for this layer, together with the opticalTable
%   field -- see loadGaAsNPOptical.m for the format).
%
% absorberIdx: index of the layer that is the deltaEc/deltaEv=0 reference (usually
%   the main absorber).
%
% junctions: struct array, length = numel(layers)-1, each containing Sn_int,
%   Sp_int (OPTIONAL, default 0/0 = ideal interface) -- defines the interface between
%   layer k and k+1.
%
% NOTE (Item 9, external review -- DELIBERATELY DEFERRED): dx=t/(N-1) (first node
% at the start of the layer, LAST node EXACTLY at the end of the layer) was tried first, but in N-layer
% chaining, since xOffset advances by EXACTLY t after each layer,
% the LAST node of layer k and the FIRST node of layer k+1 coincide at the SAME x
% -> h=0 -> eps/(h*hAvg)=Inf -> Newton/Thomas breaks down (see Task #7 debug).
% To preserve EXACTLY the same structure (and regression) as common/buildDevice.m,
% dx = t/N is used here as well (the same "one cell short" behavior as the old code,
% which also naturally preserves the NONZERO iE/iA gap between ETL/Abs that
% interfaceRecombination.m needs). The proper fix of Item 9
% (N-layer shared boundary node design: which layer's
% properties it takes, how R_int is split) was DEFERRED as a separate task
% that needs careful design.

c = constants();
NL = numel(layers);
if nargin < 3 || isempty(junctions)
  junctions = repmat(struct('Sn_int',0,'Sp_int',0), 1, NL-1);
end
if nargin < 4 || isempty(opts), opts = struct(); end
if ~isfield(opts,'Sback') || isempty(opts.Sback), opts.Sback = 1e4; end

defL = struct('Nd',0,'Na',0,'Brad',0,'CnAuger',0,'CpAuger',0, ...
  'alphaMax',0,'edgeWidth',0.05, ...
  'useOpticalTable',false,'opticalTable',[],'Nc',[],'Nv',[], ...
  'srhCenters',[]);
for k = 1:NL
  fn = fieldnames(defL);
  for i = 1:numel(fn)
    if ~isfield(layers(k), fn{i}) || isempty(layers(k).(fn{i}))
      layers(k).(fn{i}) = defL.(fn{i});
    end
  end
  if isempty(layers(k).Nc), layers(k).Nc = layers(k).NcNv; end
  if isempty(layers(k).Nv), layers(k).Nv = layers(k).NcNv; end
end

Ntotal = sum([layers.N]);
x=zeros(Ntotal,1); eps_=zeros(Ntotal,1); Eg=zeros(Ntotal,1);
Nd=zeros(Ntotal,1); Na=zeros(Ntotal,1); muN=zeros(Ntotal,1); muP=zeros(Ntotal,1);
tau=zeros(Ntotal,1); ni=zeros(Ntotal,1); chi=zeros(Ntotal,1);
Nc=zeros(Ntotal,1); Nv=zeros(Ntotal,1);
Brad=zeros(Ntotal,1); CnAuger=zeros(Ntotal,1); CpAuger=zeros(Ntotal,1);
alphaMaxArr=zeros(Ntotal,1); edgeWidthArr=zeros(Ntotal,1);
layerIdxOfNode=zeros(Ntotal,1); useTableArr=false(Ntotal,1);
% Per-NODE optical table (not per-device): previously only a single
% device-level dev.alphaTable existed (set separately by caller scripts,
% e.g. dev.alphaTable=opt), and singleWavelengthGenML.m applied THAT ONE
% table to every useOpticalTable=true node regardless of which layer it
% belonged to -- silently wrong for any device mixing two different real
% materials (e.g. a GaP window + a GaAsNP absorber both on real optical
% data), confirmed this session (window ended up absorbing like the
% absorber's own material, not like real GaP). optTableArr stores each
% node's OWN layer's opticalTable struct.
optTableArr = cell(Ntotal,1);

% Determine the global number of explicit SRH centers. Layers without an
% explicit definition receive the legacy single mid-gap center. A uniform
% matrix width keeps the continuity equations vectorized.
nCenters=1;
for k=1:NL
  if ~isempty(layers(k).srhCenters)
    nCenters=max(nCenters,numel(layers(k).srhCenters));
  end
end
srhTauN=Inf(Ntotal,nCenters); srhTauP=Inf(Ntotal,nCenters);
srhN1=zeros(Ntotal,nCenters); srhP1=zeros(Ntotal,nCenters);

xOffset = 0; node = 0;
junctionNodes = zeros(1,NL-1); % LAST node of layer k (interface k, between k and k+1)
layerStartNode = zeros(1,NL);
layerEndNode = zeros(1,NL);
for k = 1:NL
  Lk = layers(k);
  layerStartNode(k) = node+1;
  if Lk.N < 2
    error('buildDeviceML: layer %d must contain at least 2 nodes', k);
  end
  dx = Lk.t / Lk.N; % same as common/buildDevice.m (see note above)
  for j = 0:(Lk.N-1)
    node = node+1;
    x(node) = xOffset + j*dx;
    eps_(node)=Lk.eps; Eg(node)=Lk.Eg; chi(node)=Lk.chi;
    Nc(node)=Lk.Nc; Nv(node)=Lk.Nv;
    Nd(node)=Lk.Nd; Na(node)=Lk.Na;
    muN(node)=Lk.muN; muP(node)=Lk.muP;
    Brad(node)=Lk.Brad;
    CnAuger(node)=Lk.CnAuger; CpAuger(node)=Lk.CpAuger;
    alphaMaxArr(node)=Lk.alphaMax; edgeWidthArr(node)=Lk.edgeWidth;
    useTableArr(node)=Lk.useOpticalTable;
    optTableArr{node}=Lk.opticalTable;
    layerIdxOfNode(node) = k;
    ni(node) = sqrt(Lk.Nc*Lk.Nv) * exp(-Eg(node) / (2*c.VT));
    tau(node) = Lk.tau;
    if isempty(Lk.srhCenters)
      srhTauN(node,1)=Lk.tau; srhTauP(node,1)=Lk.tau;
      srhN1(node,1)=ni(node); srhP1(node,1)=ni(node);
    else
      for jc=1:numel(Lk.srhCenters)
        ctr=Lk.srhCenters(jc);
        if ~isfield(ctr,'tauN') || ~isfield(ctr,'tauP') || ...
           ~isfield(ctr,'EtRelEi_eV')
          error('buildDeviceML:SRHCenterFields', ...
            'Each SRH center requires tauN, tauP, and EtRelEi_eV.');
        end
        srhTauN(node,jc)=ctr.tauN; srhTauP(node,jc)=ctr.tauP;
        srhN1(node,jc)=ni(node)*exp(ctr.EtRelEi_eV/c.VT);
        srhP1(node,jc)=ni(node)*exp(-ctr.EtRelEi_eV/c.VT);
      end
    end
  end
  layerEndNode(k) = node;
  xOffset = xOffset + Lk.t;
  if k < NL
    junctionNodes(k) = node;
  end
end

h = diff(x);

chiAbs = layers(absorberIdx).chi;
EgAbs = layers(absorberIdx).Eg;
NcAbs = layers(absorberIdx).Nc;
NvAbs = layers(absorberIdx).Nv;
chiRel = chi - chiAbs;
psiOffsetN = chiRel + c.VT*log(Nc/NcAbs);
psiOffsetP = chiRel + (Eg-EgAbs) - c.VT*log(Nv/NvAbs);
% Compatibility aliases for older diagnostic scripts. The transport code
% uses psi+psiOffset; old code used psi-deltaE.
deltaEc = -psiOffsetN;
deltaEv = -psiOffsetP;

% isAbs: kept for compatibility with old 2-region code (solveElectronContinuity
% etc. do NOT use this field DIRECTLY -- only singleWavelengthGenML.m uses it).
% Here ALL nodes are considered potentially absorbing (if alphaMax=0 or
% useOpticalTable=false the contribution is zero anyway); for multilayer optics
% isAbs was made independent of a single "which region" distinction.
isAbs = double(useTableArr | (alphaMaxArr > 0));

Sn_int_arr = arrayfun(@(j) j.Sn_int, junctions);
Sp_int_arr = arrayfun(@(j) j.Sp_int, junctions);

% If there is a single interface (2 layers), we also fill the BACKWARD-COMPATIBLE
% single-junction fields (N_ETL) used by common/interfaceRecombination.m --
% so that common/ files work UNCHANGED in the 2-layer special case
% without modification (see Task #7 regression validation).
dev = struct('N',Ntotal,'x',x,'h',h,'eps',eps_,'Eg',Eg,'Nd',Nd,'Na',Na, ...
  'muN',muN,'muP',muP,'tau',tau,'isAbs',isAbs,'ni',ni,'chi',chi, ...
  'Nc',Nc,'Nv',Nv,'chiRef',chiAbs,'EgRef',EgAbs, ...
  'psiOffsetN',psiOffsetN,'psiOffsetP',psiOffsetP, ...
  'deltaEc',deltaEc,'deltaEv',deltaEv,'Brad',Brad, ...
  'CnAuger',CnAuger,'CpAuger',CpAuger, ...
  'srhTauN',srhTauN,'srhTauP',srhTauP,'srhN1',srhN1,'srhP1',srhP1, ...
  'alphaMaxNode',alphaMaxArr,'edgeWidthNode',edgeWidthArr,'useOpticalTable',useTableArr, ...
  'opticalTableNode',{optTableArr}, ...
  'layerIdxOfNode',layerIdxOfNode,'layerStartNode',layerStartNode,'layerEndNode',layerEndNode, ...
  'junctionNodes',junctionNodes,'Sn_int_j',Sn_int_arr,'Sp_int_j',Sp_int_arr, ...
  'absorberIdx',absorberIdx,'Sback',opts.Sback);

if NL >= 2
  dev.N_ETL = junctionNodes(1); % backward compatibility with common/interfaceRecombination.m (first interface)
  dev.Sn_int = Sn_int_arr(1);
  dev.Sp_int = Sp_int_arr(1);
else
  dev.N_ETL = 0; dev.Sn_int = 0; dev.Sp_int = 0;
end
end
