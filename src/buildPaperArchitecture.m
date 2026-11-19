function [dev, meta] = buildPaperArchitecture(architectureCode, opts)
%BUILDPAPERARCHITECTURE Canonical builder for the seven paper architectures.
%
% [dev,meta] = buildPaperArchitecture(code) builds one controlled-comparison
% device defined by docs/ARCHITECTURE_CONFIGURATION_MATRIX.md. Supported
% codes are HJ-0, HJ-U, HJ-G, F-HJ, RBSF-GaP, RBSF-GaNP, and RHJ.
% Name-value options include solveEquilibrium, frontIRV, rearIRV, and Sback.
%
% The builder deliberately centralizes every common parameter. Architecture
% branches may change only layer sequence, doping profile, and rear material.
% Existing numbered case scripts are legacy and are not called here.

arguments
  architectureCode (1,1) string
  opts.solveEquilibrium (1,1) logical = true
  opts.frontIRV (1,1) double {mustBeNonnegative} = 8910
  opts.rearIRV (1,1) double {mustBeNonnegative} = 8910
  opts.Sback (1,1) double {mustBeNonnegative} = 1e5
  opts.BradAbs (1,1) double {mustBeNonnegative} = 0
  opts.BradGaP (1,1) double {mustBeNonnegative} = 0
  opts.BradGaNP (1,1) double {mustBeNonnegative} = 0
  opts.CnAugerAbs (1,1) double {mustBeNonnegative} = 0
  opts.CpAugerAbs (1,1) double {mustBeNonnegative} = 0
  opts.CnAugerGaP (1,1) double {mustBeNonnegative} = 0
  opts.CpAugerGaP (1,1) double {mustBeNonnegative} = 0
  opts.CnAugerGaNP (1,1) double {mustBeNonnegative} = 0
  opts.CpAugerGaNP (1,1) double {mustBeNonnegative} = 0
  opts.meshScale (1,1) double {mustBePositive} = 1
  opts.GaPChi_eV (1,1) double {mustBePositive} = 3.8
  opts.GaPTau_s (1,1) double {mustBePositive} = 1e-9
  opts.AbsorberDoping (1,1) double {mustBePositive} = 1e17   % GaAsPN absorber/base doping (p for HJ/F-HJ/RBSF, n for RHJ); default frozen at 1e17 cm-3
  opts.AbsorberTau_s (1,1) double {mustBePositive} = 0.6e-6   % GaAsPN SRH lifetime; default frozen at the Ilahi-derived nominal
  opts.AbsorberMuN (1,1) double {mustBePositive} = 500        % GaAsPN electron mobility (assumed nominal), cm2/Vs
  opts.HomoRearContactNa (1,1) double {mustBePositive} = 1e19
  opts.HomoRearContactNodes (1,1) double {mustBeInteger,mustBeGreaterThanOrEqual(opts.HomoRearContactNodes,2)} = 20
  opts.AbsorberOpticalShift_eV (1,1) double = 0     % rigid shift of the GaAsPN optical table in energy (referee check; 0 = frozen production)
  opts.AbsorberSRHCenters = []                      % struct array (tauN,tauP,EtRelEi_eV) of SRH centres of the GaAsPN absorber layers (referee check, MD421-calibrated regime); [] = single centre with opts.AbsorberTau_s (frozen production)
  opts.AbsorberEgShift_eV (1,1) double = 0          % shift of the GaAsPN electrical bandgap (referee check; use together with AbsorberOpticalShift_eV and a consistent BradAbs; 0 = frozen production)
  opts.AbsorberEgShiftInCB (1,1) logical = true      % true: the Eg shift is taken in the conduction band (chi_abs = 4.0 - shift, valence band fixed); false: chi fixed (valence band moves)
  opts.AbsorberThicknessScale (1,1) double {mustBePositive} = 1   % scales every GaAsPN absorber/base thickness and its node count (referee check; 1 = frozen production)
  opts.PinBackground (1,1) double = 1e15   % PIN only: background doping of the unintentionally doped GaAsPN i-layer, cm-3 (>0 n-type, <0 p-type)
  opts.InterfaceSlabNodes (1,1) double {mustBeInteger,mustBeNonnegative} = 0   % referee check R1-2: >0 splits a slab of InterfaceSlab_nm off each side of every heterointerface and meshes it with this many nodes (0 = frozen production mesh)
  opts.InterfaceSlab_nm (1,1) double {mustBePositive} = 5
  opts.AbsorberEps (1,1) double {mustBePositive} = 12.9   % GaAsPN relative permittivity (round-5 review: Vegard ~11.1; 12.9 = frozen production)
  opts.AbsorberNc (1,1) double {mustBePositive} = 1e19    % GaAsPN effective densities of states (assumption; frozen production 1e19/1e19)
  opts.AbsorberNv (1,1) double {mustBePositive} = 1e19
  opts.BGN_eV (1,1) double {mustBeNonnegative} = 0        % round-6 check: band-gap narrowing of heavily doped layers (0 = frozen production)
  opts.BGNThreshold (1,1) double {mustBePositive} = 1e18  % doping above which BGN_eV is applied
  opts.IntRecLowGap (1,1) logical = false                 % round-6 check: interface recombination with lower-gap-side densities (false = geometric mean, production)
  opts.AbsorberNFraction (1,1) double {mustBePositive} = 0.03   % nitrogen fraction of the Si-matched GaAsPN absorber (composition check; 0.03 = frozen production). Eg from BAC; the gap change relative to 3% N is taken in the conduction band (chi_abs = 4.0 - dEg)
  opts.AbsorberOpticalN (1,1) double {mustBePositive} = 3.0     % percent N of the Zou-Goodnick n,k table (2.5,3,...,5; 3.0 = frozen production)
end

code = canonicalCode(architectureCode);

gap = gapParams();
absBand = gaasnpBandgap(opts.AbsorberNFraction);
dEgN = absBand.Eg - gaasnpBandgap(0.03).Eg;   % 0 for the production composition
gaPOpt = loadGaPOptical();
absOpt = loadGaAsNPOptical(opts.AbsorberOpticalN);
if opts.AbsorberOpticalShift_eV ~= 0
  absOpt.E = absOpt.E + opts.AbsorberOpticalShift_eV;
  absOpt.Emin = min(absOpt.E); absOpt.Emax = max(absOpt.E);
end
tS = opts.AbsorberThicknessScale;
gaNP = ganp21Params();

% Shared material definitions. Values and provenance are frozen in the
% architecture matrix; uncertain quantities must be varied through a later
% paired uncertainty runner, not tuned separately inside this function.
matGaP = struct('Eg',2.26,'chi',opts.GaPChi_eV,'eps',11.1,'Nc',1.8e19,'Nv',1.9e19, ...
  'muN',100,'muP',80,'tau',opts.GaPTau_s,'Brad',opts.BradGaP, ...
  'CnAuger',opts.CnAugerGaP,'CpAuger',opts.CpAugerGaP);
chiAbs = 4.0 - dEgN; if opts.AbsorberEgShiftInCB, chiAbs = chiAbs - opts.AbsorberEgShift_eV; end
matAbs = struct('Eg',absBand.Eg + opts.AbsorberEgShift_eV,'chi',chiAbs,'eps',opts.AbsorberEps,'Nc',opts.AbsorberNc,'Nv',opts.AbsorberNv, ...
  'muN',opts.AbsorberMuN,'muP',112,'tau',opts.AbsorberTau_s,'Brad',opts.BradAbs, ...
  'CnAuger',opts.CnAugerAbs,'CpAuger',opts.CpAugerAbs,'srhCenters',opts.AbsorberSRHCenters);
% GaNP transport values are explicit assumptions inherited from the GaP
% rear layer; meta.provenance records that they are not GaNP measurements.
matGaNP = struct('Eg',gaNP.Eg,'chi',4.10,'eps',11.1,'Nc',matGaP.Nc,'Nv',matGaP.Nv, ...
  'muN',matGaP.muN,'muP',matGaP.muP,'tau',1e-9,'Brad',opts.BradGaNP, ...
  'CnAuger',opts.CnAugerGaNP,'CpAuger',opts.CpAugerGaNP);

layers = repmat(emptyLayer(),1,0);
layerNames = strings(1,0);
junctions = repmat(struct('Sn_int',0,'Sp_int',0),1,0);

% Every architecture starts with the same optically active GaP window and
% the same physical GaP/GaAsPN interface recombination velocity.
if code == "RBSF-GaP-n"
  % Polarity check only (not a paper architecture): p+ window for the n-type front-junction mirror of RBSF-GaP.
  [layers,layerNames] = appendLayer(layers,layerNames, ...
    makeLayer("p+-GaP window","GaP",30e-7,30,matGaP,0,1e19,true,gaPOpt));
else
  [layers,layerNames] = appendLayer(layers,layerNames, ...
    makeLayer("n+-GaP window","GaP",30e-7,30,matGaP,1e19,0,true,gaPOpt));
end

switch code
  case {"HJ-0","HJ-U","HJ-G"}
    [layers,layerNames] = appendLayer(layers,layerNames, ...
      makeLayer("n-GaAsPN emitter","GaAsPN",150e-7,60,matAbs,1e18,0,true,absOpt));
    junctions(end+1) = makeJunction(opts.frontIRV);

    if code == "HJ-0"
      % True BSF-free control: the uniformly doped base reaches the rear
      % surface directly. Its extra 100 nm preserves the common 1500 nm
      % GaAsPN optical-thickness budget after removing the p+ contact layer.
      nSlices = 27;
      baseThickness = 1350e-7*tS;
    else
      nSlices = 25;
      baseThickness = 1250e-7*tS;
    end
    sliceThickness = baseThickness/nSlices;
    for k = 1:nSlices
      xMid = (k-0.5)/nSlices;
      if code == "HJ-G"
        Na = opts.AbsorberDoping*exp(log(10)*xMid);
      else
        Na = opts.AbsorberDoping;
      end
      [layers,layerNames] = appendLayer(layers,layerNames, ...
        makeLayer(sprintf("p-GaAsPN base %02d",k),"GaAsPN", ...
          sliceThickness,8,matAbs,0,Na,true,absOpt));
      % The n/p transition is a homojunction and numerical base slices are
      % not physical recombination interfaces.
      junctions(end+1) = makeJunction(0);
    end
    if code ~= "HJ-0"
      [layers,layerNames] = appendLayer(layers,layerNames, ...
        makeLayer("p+-GaAsPN contact","GaAsPN",100e-7,opts.HomoRearContactNodes,matAbs,0,opts.HomoRearContactNa,true,absOpt));
      junctions(end+1) = makeJunction(0);
    end

  case "F-HJ"
    [layers,layerNames] = appendLayer(layers,layerNames, ...
      makeLayer("p-GaAsPN absorber","GaAsPN",1400e-7*tS,max(2,round(224*tS)),matAbs,0,opts.AbsorberDoping,true,absOpt));
    junctions(end+1) = makeJunction(opts.frontIRV);
    [layers,layerNames] = appendLayer(layers,layerNames, ...
      makeLayer("p+-GaAsPN contact","GaAsPN",100e-7,opts.HomoRearContactNodes,matAbs,0,opts.HomoRearContactNa,true,absOpt));
    junctions(end+1) = makeJunction(0);

  case "RBSF-GaP"
    [layers,layerNames] = appendLayer(layers,layerNames, ...
      makeLayer("p-GaAsPN absorber","GaAsPN",1500e-7*tS,max(2,round(240*tS)),matAbs,0,opts.AbsorberDoping,true,absOpt));
    junctions(end+1) = makeJunction(opts.frontIRV);
    [layers,layerNames] = appendLayer(layers,layerNames, ...
      makeLayer("p+-GaP rear BSF","GaP",100e-7,20,matGaP,0,1e19,false,[]));
    junctions(end+1) = makeJunction(opts.rearIRV);

  case "RBSF-GaNP"
    [layers,layerNames] = appendLayer(layers,layerNames, ...
      makeLayer("p-GaAsPN absorber","GaAsPN",1500e-7*tS,max(2,round(240*tS)),matAbs,0,opts.AbsorberDoping,true,absOpt));
    junctions(end+1) = makeJunction(opts.frontIRV);
    [layers,layerNames] = appendLayer(layers,layerNames, ...
      makeLayer("p+-GaNP rear BSF","GaNP",100e-7,20,matGaNP,0,1e19,false,[]));
    junctions(end+1) = makeJunction(opts.rearIRV);

  case "PIN"
    % n+-GaP / i-GaAsPN / p+-GaP double heterostructure (the structure of the
    % grown GaAsPN cells): the RBSF-GaP stack with the absorber doping replaced
    % by a low background, so that the field of the n+/p+ pair extends over the
    % i-layer (referee check R2-4 / R5-M2).
    NdI = max(opts.PinBackground,0); NaI = max(-opts.PinBackground,0);
    [layers,layerNames] = appendLayer(layers,layerNames, ...
      makeLayer("i-GaAsPN absorber","GaAsPN",1500e-7*tS,max(2,round(240*tS)),matAbs,NdI,NaI,true,absOpt));
    junctions(end+1) = makeJunction(opts.frontIRV);
    [layers,layerNames] = appendLayer(layers,layerNames, ...
      makeLayer("p+-GaP rear BSF","GaP",100e-7,20,matGaP,0,1e19,false,[]));
    junctions(end+1) = makeJunction(opts.rearIRV);

  case "RBSF-GaP-n"
    % n-type absorber with the junction at the FRONT (p+-GaP / n-GaAsPN / n+-GaP rear BSF): separates the absorber type
    % from the junction position in the RBSF-GaP versus RHJ comparison (check, 2026-11-19).
    [layers,layerNames] = appendLayer(layers,layerNames, ...
      makeLayer("n-GaAsPN absorber","GaAsPN",1500e-7*tS,max(2,round(240*tS)),matAbs,opts.AbsorberDoping,0,true,absOpt));
    junctions(end+1) = makeJunction(opts.frontIRV);
    [layers,layerNames] = appendLayer(layers,layerNames, ...
      makeLayer("n+-GaP rear BSF","GaP",100e-7,20,matGaP,1e19,0,false,[]));
    junctions(end+1) = makeJunction(opts.rearIRV);

  case "RHJ"
    [layers,layerNames] = appendLayer(layers,layerNames, ...
      makeLayer("n-GaAsPN absorber","GaAsPN",1500e-7*tS,max(2,round(240*tS)),matAbs,opts.AbsorberDoping,0,true,absOpt));
    junctions(end+1) = makeJunction(opts.frontIRV);
    [layers,layerNames] = appendLayer(layers,layerNames, ...
      makeLayer("p+-GaP rear junction","GaP",100e-7,20,matGaP,0,1e19,false,[]));
    junctions(end+1) = makeJunction(opts.rearIRV);
end

% Mesh scaling is applied after the physical architecture has been assembled,
% so thicknesses, compositions, doping profiles, and interface locations are
% identical in every discretization-convergence run.
for k = 1:numel(layers)
  layers(k).N = max(2,round(layers(k).N*opts.meshScale));
end

absorberIdx = 2;
if opts.InterfaceSlabNodes > 0
  % Local refinement at heterointerfaces: the band step sits in the segment between the last node of the left layer and the first
  % node of the right layer, whose width is the node spacing of the left layer; a thin slab of the same material with fine nodes on
  % each side shrinks that segment (5 nm / 25 nodes -> 0.2 nm). The slabs are numerical sub-layers (no interface recombination
  % between a slab and the rest of its layer); the physical interfaces keep their recombination velocity.
  ts = opts.InterfaceSlab_nm*1e-7; ns = opts.InterfaceSlabNodes;
  nL = numel(layers); newL = repmat(emptyLayer(),1,0); newNames = strings(1,0); newJ = repmat(struct('Sn_int',0,'Sp_int',0),1,0);
  for k = 1:nL
    L = layers(k);
    leftH = k > 1 && layers(k-1).material ~= L.material;
    rightH = k < nL && layers(k+1).material ~= L.material;
    tRem = L.t - ts*(leftH + rightH);
    if tRem <= ts, leftH = false; rightH = false; tRem = L.t; end   % too thin to split
    pieces = repmat(emptyLayer(),1,0);
    if leftH, S = L; S.t = ts; S.N = ns; S.name = L.name + " [slab]"; pieces(end+1) = S; end %#ok<AGROW>
    M = L; M.t = tRem; M.N = max(2,round(L.N*tRem/L.t)); pieces(end+1) = M;
    if rightH, S = L; S.t = ts; S.N = ns; S.name = L.name + " [slab]"; pieces(end+1) = S; end %#ok<AGROW>
    if k == 2, absorberIdx = numel(newL) + 1 + leftH; end
    for q = 1:numel(pieces)
      newL(end+1) = pieces(q); newNames(end+1) = pieces(q).name; %#ok<AGROW>
      if q < numel(pieces), newJ(end+1) = makeJunction(0); %#ok<AGROW>
      elseif k < nL, newJ(end+1) = junctions(k); end %#ok<AGROW>
    end
  end
  layers = newL; layerNames = newNames; junctions = newJ;
end
if opts.BGN_eV > 0
  % round-6 review check: rigid band-gap narrowing of every layer doped >= opts.BGNThreshold; p-type: valence band up (chi fixed),
  % n-type: conduction band down (chi + dEg). Optical tables, B_rad and lifetimes unchanged.
  for k = 1:numel(layers)
    Ndop = max(layers(k).Nd, layers(k).Na);
    if Ndop >= opts.BGNThreshold
      layers(k).Eg = layers(k).Eg - opts.BGN_eV;
      if layers(k).Nd > layers(k).Na, layers(k).chi = layers(k).chi + opts.BGN_eV; end
    end
  end
end
if opts.IntRecLowGap, intRecMode = "lowgap"; else, intRecMode = "geomean"; end
dev = buildDeviceML(layers,absorberIdx,junctions,struct('Sback',opts.Sback));
dev.intRecMode = intRecMode;
if opts.solveEquilibrium
  dev.psiEq = solveEquilibrium(dev);
else
  dev.psiEq = [];
end

isGaAsPN = arrayfun(@(x) string(x.material)=="GaAsPN",layers);
meta = struct();
meta.schemaVersion = "paper-architecture-v1";
meta.architectureCode = code;
meta.architectureName = architectureName(code);
meta.legacyAlias = legacyAlias(code);
meta.layerNames = layerNames;
meta.layers = layers;
meta.junctions = junctions;
meta.totalGaAsPNThickness_cm = sum([layers(isGaAsPN).t]);
meta.totalGaAsPNThickness_nm = meta.totalGaAsPNThickness_cm*1e7;
meta.frontIRV_cmps = opts.frontIRV;
meta.rearIRV_cmps = opts.rearIRV;
meta.Sback_cmps = opts.Sback;
meta.meshScale = opts.meshScale;
meta.totalMeshNodes = sum([layers.N]);
meta.GaPChi_eV = opts.GaPChi_eV;
meta.GaPTau_s = opts.GaPTau_s;
meta.HomoRearContactNa_cm3 = opts.HomoRearContactNa;
meta.HomoRearContactNodes = opts.HomoRearContactNodes;
meta.recombination = struct('BradAbs_cm3ps',opts.BradAbs, ...
  'BradGaP_cm3ps',opts.BradGaP,'BradGaNP_cm3ps',opts.BradGaNP, ...
  'CnAugerAbs_cm6ps',opts.CnAugerAbs,'CpAugerAbs_cm6ps',opts.CpAugerAbs, ...
  'CnAugerGaP_cm6ps',opts.CnAugerGaP,'CpAugerGaP_cm6ps',opts.CpAugerGaP, ...
  'CnAugerGaNP_cm6ps',opts.CnAugerGaNP,'CpAugerGaNP_cm6ps',opts.CpAugerGaNP);
meta.absorberComposition = struct('N_fraction',absBand.N_fraction, ...
  'P_fraction',absBand.P_fraction,'As_fraction',absBand.As_fraction, ...
  'Eg_eV',absBand.Eg);
meta.gaNPComposition = gaNP;
meta.provenance = struct( ...
  'GaP_transport',"literature-guided nominal", ...
  'GaAsPN_muN',"assumed; uncertainty required", ...
  'GaAsPN_muP_tau',"Ilahi-derived nominal", ...
  'GaAsPN_NcNv_chi',"assumed; uncertainty required", ...
  'GaNP_transport',"GaP-borrowed assumption; uncertainty required");

dev.paperArchitecture = rmfield(meta,{'layers','junctions'});
end

function code = canonicalCode(raw)
key = upper(strrep(strtrim(raw),"_","-"));
switch key
  case {"HJ-0","HJ0","NO-BSF"}, code = "HJ-0";
  case {"HJ-U","D1","#49-D1","49-D1"}, code = "HJ-U";
  case {"HJ-G","D2","#49-D2","49-D2"}, code = "HJ-G";
  case {"F-HJ","#50","50"}, code = "F-HJ";
  case {"RBSF-GAP","#57","57"}, code = "RBSF-GaP";
  case {"RBSF-GANP","#58","58"}, code = "RBSF-GaNP";
  case {"RHJ"}, code = "RHJ";
  case {"RBSF-GAP-N"}, code = "RBSF-GaP-n";
  case {"PIN","P-I-N","NIP"}, code = "PIN";
  otherwise
    error('buildPaperArchitecture:UnknownArchitecture', ...
      'Unknown architecture "%s". Use HJ-0, HJ-U, HJ-G, F-HJ, RBSF-GaP, RBSF-GaNP, or RHJ.',raw);
end
end

function name = architectureName(code)
switch code
  case "HJ-0", name = "BSF-free homojunction control";
  case "HJ-U", name = "Uniform-doped homojunction";
  case "HJ-G", name = "Graded-doping homojunction";
  case "F-HJ", name = "Front heterojunction";
  case "RBSF-GaP", name = "GaP rear-BSF heterostructure";
  case "RBSF-GaNP", name = "GaNP rear-BSF heterostructure";
  case "RHJ", name = "Rear heterojunction";
  case "PIN", name = "n+-GaP/i-GaAsPN/p+-GaP double heterostructure";
  case "RBSF-GaP-n", name = "n-type front-junction mirror of RBSF-GaP (p+-GaP/n-GaAsPN/n+-GaP)";
end
end

function alias = legacyAlias(code)
switch code
  case "HJ-0", alias = "new control";
  case "HJ-U", alias = "#49-D1";
  case "HJ-G", alias = "#49-D2";
  case "F-HJ", alias = "#50";
  case "RBSF-GaP", alias = "#57";
  case "RBSF-GaNP", alias = "#58";
  case "RHJ", alias = "RHJ";
  case "PIN", alias = "p-i-n";
  case "RBSF-GaP-n", alias = "n-mirror";
end
end

function J = makeJunction(irv)
J = struct('Sn_int',irv,'Sp_int',irv);
end

function L = emptyLayer()
L = struct('name',"",'material',"",'t',0,'N',0,'Eg',0,'chi',0,'eps',0, ...
  'Nd',0,'Na',0,'muN',0,'muP',0,'tau',0,'Nc',0,'Nv',0,'NcNv',0, ...
  'Brad',0,'CnAuger',0,'CpAuger',0,'alphaMax',0,'edgeWidth',0.05,'useOpticalTable',false, ...
  'opticalTable',[],'srhCenters',[]);
end

function L = makeLayer(name,material,t,N,mat,Nd,Na,useOpticalTable,opticalTable)
L = emptyLayer();
L.name = string(name); L.material = string(material);
L.t=t; L.N=N; L.Eg=mat.Eg; L.chi=mat.chi; L.eps=mat.eps;
L.Nd=Nd; L.Na=Na; L.muN=mat.muN; L.muP=mat.muP; L.tau=mat.tau;
L.Nc=mat.Nc; L.Nv=mat.Nv; L.NcNv=sqrt(mat.Nc*mat.Nv); L.Brad=mat.Brad;
L.CnAuger=mat.CnAuger; L.CpAuger=mat.CpAuger;
if isfield(mat,'srhCenters') && ~isempty(mat.srhCenters), L.srhCenters = mat.srhCenters; end
L.useOpticalTable=useOpticalTable; L.opticalTable=opticalTable;
end

function [layers,names] = appendLayer(layers,names,L)
layers(end+1)=L;
names(end+1)=L.name;
end
