function dev = buildDeviceML(layers, absorberIdx, junctions, opts)
% Genel N-katmanli cihaz kurucu. buildDevice.m'nin (common/, 2-bolgeli
% ETL/Abs) genellestirilmis hali -- var olan dosyaya DOKUNMADAN, ayri bir
% modul olarak. Ciktisi (dev struct) common/'daki solveElectronContinuity.m
% ve solveHoleContinuity.m tarafindan DOGRUDAN kullanilabilir (ayni alan
% adlari/formatlari), cunku bu iki dosya zaten herhangi bir bolge sayisina
% genellenebilir sekilde (sadece per-node dizilerle) yazilmis.
%
% layers: struct dizisi, her biri:
%   t (kalinlik, cm), N (izgara nokta sayisi), Eg, chi, eps, Nd, Na
%   (katkilama, birindeki OPSIYONEL, varsayilan 0), muN, muP, tau, NcNv,
%   Brad (OPSIYONEL, varsayilan 0), alphaMax, edgeWidth (OPSIYONEL,
%   lojistik-kenar sogurma modeli icin; kullanilmazsa 0/0.05 varsayilan --
%   yani o katman optik olarak SEFFAF), useOpticalTable (OPSIYONEL, true
%   ise bu katman icin gercek n(E),k(E) tablosu kullanilir, opticalTable
%   alaniyla birlikte -- bkz. loadGaAsNPOptical.m formatinda).
%
% absorberIdx: deltaEc/deltaEv=0 referansi olan katmanin indeksi (genelde
%   ana absorber).
%
% junctions: struct dizisi, uzunluk = numel(layers)-1, her biri Sn_int,
%   Sp_int (OPSIYONEL, varsayilan 0/0 = ideal arayuz) icerir -- katman k
%   ile k+1 arasindaki arayuzu tanimlar.
%
% NOT (Madde 9, dis inceleme -- BILEREK ERTELENDI): dx=t/(N-1) (ilk dugum
% katman basinda, SON dugum TAM katman sonunda) once denendi, ama N-katmanli
% zincirlemede xOffset her katmandan sonra TAM t kadar ilerledigi icin
% katman k'nin SON dugumu ile katman k+1'in ILK dugumu AYNI x'te cakisiyor
% -> h=0 -> eps/(h*hAvg)=Inf -> Newton/Thomas coker (bkz. Gorev #7 debug).
% common/buildDevice.m ile BIREBIR ayni yapiyi (ve regresyonu) korumak icin
% burada da dx = t/N kullaniliyor (eski kod ile ayni "bir hucre kisa kalma"
% davranisi, ETL/Abs arasinda interfaceRecombination.m'nin ihtiyac duydugu
% SIFIR-OLMAYAN iE/iA araligini da dogal olarak koruyor). Madde 9'un asil
% duzeltmesi (N-katmanli ortak sinir dugumu tasarimi: hangi katmanin
% ozelliklerini alacagi, R_int'in nasil bolunecegi) AYRI, dikkatli
% tasarlanmasi gereken bir gorev olarak ERTELENDI.

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
junctionNodes = zeros(1,NL-1); % katman k'nin SON dugumu (arayuz k, k ile k+1 arasinda)
layerStartNode = zeros(1,NL);
layerEndNode = zeros(1,NL);
for k = 1:NL
  Lk = layers(k);
  layerStartNode(k) = node+1;
  if Lk.N < 2
    error('buildDeviceML: katman %d en az 2 dugum icermeli', k);
  end
  dx = Lk.t / Lk.N; % common/buildDevice.m ile ayni (bkz. yukaridaki not)
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

% isAbs: eski 2-bolgeli kodla uyumluluk icin tutuluyor (solveElectronContinuity
% vb. bu alani DOGRUDAN kullanmiyor -- sadece singleWavelengthGenML.m kullanir).
% Burada TUM dugumler potansiyel olarak soguruyor sayilir (alphaMax=0 veya
% useOpticalTable=false ise zaten katki sifir olur), coklu-katman optigi
% icin isAbs tek bir "hangi bolge" ayrimindan bagimsizlastirildi.
isAbs = double(useTableArr | (alphaMaxArr > 0));

Sn_int_arr = arrayfun(@(j) j.Sn_int, junctions);
Sp_int_arr = arrayfun(@(j) j.Sp_int, junctions);

% Tek bir arayuz varsa (2 katman), common/interfaceRecombination.m ile
% GERIYE-UYUMLU tek-junction alanlarini da (N_ETL) dolduruyoruz --
% common/ dosyalari degistirilmeden 2-katmanli ozel durumda AYNEN
% calismasi icin (bkz. Gorev #7 regresyon dogrulamasi).
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
  dev.N_ETL = junctionNodes(1); % common/interfaceRecombination.m ile geriye-uyumluluk (ilk arayuz)
  dev.Sn_int = Sn_int_arr(1);
  dev.Sp_int = Sp_int_arr(1);
else
  dev.N_ETL = 0; dev.Sn_int = 0; dev.Sp_int = 0;
end
end
