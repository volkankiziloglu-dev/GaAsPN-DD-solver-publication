function results = run_geisz2002_staged_calibration()
%RUN_GEISZ2002_STAGED_CALIBRATION Staged MD421 calibration.
%
% Separates:
%   1) an effective optical/collection scale constrained by Jsc and IQE,
%   2) two-center SRH candidates constrained by Voc,
%   3) terminal Rs/Rsh constrained primarily by FF.
%
% The effective collection scale is not presented as an absorption
% coefficient measurement: it represents missing reflectance, emitter
% absorption and collection losses until composition-specific optical
% stacks and one-sun IQE become available.
assertCanonicalPath();
root=fileparts(fileparts(mfilename('fullpath')));
gap=gapParams(); opt0=loadGaAsNPOptical(3.0);
am15=loadRealAM15G(); Pin=spectrumPower(am15);
target=struct('Jsc',5.8e-3,'Voc',1.12,'FF',0.60,'PCE',3.8976);

% Stage 1: baseline result gives a measured-current collection scale.
baseline=simulateMD421(opt0,[],1.0,8910,gap,am15);
collection0=min(max(target.Jsc/baseline.metrics.Jsc,0.10),1.0);
collectionGrid=unique(max(0.10,min(1.0,collection0.*[0.925 1 1.075])));

% Stage 2: physically explicit two-center candidates. Center 1 is midgap;
% center 2 is 0.25 eV above Ei with asymmetric capture times. All values
% are screening priors, not independently measured trap assignments.
tauMid=[2e-9 1e-9 5e-10];
rows=cell(0,19); curves=cell(0,1); configs=cell(0,1);
for gScale=collectionGrid
 for tm=tauMid
  ctr(1)=struct('tauN',tm,'tauP',tm,'EtRelEi_eV',0);
  ctr(2)=struct('tauN',30*tm,'tauP',3*tm,'EtRelEi_eV',0.25);
  try
   sim=simulateMD421(opt0,ctr,gScale,8910,gap,am15);
   [Rs,Rsh,mTerm,termLoss]=fitTerminal(sim.pts,Pin,target);
   m=sim.metrics;
   loss=((m.Jsc-target.Jsc)/(0.10*target.Jsc))^2 + ...
        ((m.Voc-target.Voc)/0.05)^2 + ...
        ((mTerm.FF-target.FF)/0.05)^2 + ...
        ((mTerm.PCE-target.PCE)/(0.15*target.PCE))^2;
   pass=abs(m.Jsc-target.Jsc)/target.Jsc<=0.10 && ...
        abs(m.Voc-target.Voc)<=0.05 && ...
        abs(mTerm.FF-target.FF)<=0.05 && ...
        abs(mTerm.PCE-target.PCE)/target.PCE<=0.15;
   rows(end+1,:)={gScale,tm,30*tm,3*tm,0.25,8910,Rs,Rsh, ...
    m.Jsc,m.Voc,m.FF,m.PCE,mTerm.Jsc,mTerm.Voc,mTerm.FF,mTerm.PCE, ...
    loss+termLoss,string(pass),string(sim.status.stopReason)}; %#ok<AGROW>
   curves{end+1}=sim.pts; configs{end+1}=ctr; %#ok<AGROW>
  catch ME
   rows(end+1,:)={gScale,tm,30*tm,3*tm,0.25,8910,NaN,NaN, ...
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,Inf,"false",string(ME.message)}; %#ok<AGROW>
   curves{end+1}=[]; configs{end+1}=ctr; %#ok<AGROW>
  end
 end
end
results=cell2table(rows,'VariableNames', ...
 {'CollectionScale','TauMid_s','TauEdgeN_s','TauEdgeP_s','EtEdgeRelEi_eV', ...
 'Sint_cms','Rs_ohmcm2','Rsh_ohmcm2','IntrinsicJsc_Acm2','IntrinsicVoc_V', ...
 'IntrinsicFF','IntrinsicPCE_percent','TerminalJsc_Acm2','TerminalVoc_V', ...
 'TerminalFF','TerminalPCE_percent','MetricLoss','AllCriteriaPass','Reason'});
results=sortrows(results,'MetricLoss');
results.Rank=(1:height(results))';
disp(results(1:min(10,height(results)),:));

outDir=fullfile(root,'results'); if ~exist(outDir,'dir'),mkdir(outDir);end
writetable(results,fullfile(outDir,'geisz2002_staged_calibration.csv'));
save(fullfile(outDir,'geisz2002_staged_calibration.mat'), ...
 'results','curves','configs','baseline','collectionGrid','tauMid','target');
end

function sim=simulateMD421(opt,centers,gScale,Sint,gap,am15)
NcGaP=1.8e19; NvGaP=1.9e19; NcAbs=1e19; NvAbs=1e19;
emitter=layer(100e-7,40,gap.Eg,gap.chi,gap.eps,1e18,0, ...
 gap.muN_doped,gap.muP_doped,gap.tau,NcGaP,NvGaP,false,[],[]);
base=layer(0.64e-4,150,1.88,4.0,12.9,0,1e16, ...
 500,112,0.6e-6,NcAbs,NvAbs,true,opt,centers);
bsf=layer(500e-7,70,gap.Eg,gap.chi,gap.eps,0,1e18, ...
 gap.muN_doped,gap.muP_doped,gap.tau,NcGaP,NvGaP,false,[],[]);
dev=setupDeviceML([emitter base bsf],2, ...
 struct('Sn_int',{Sint,Sint},'Sp_int',{Sint,Sint}),struct('Sback',1e5));
dev.alphaTable=opt;
G=gScale*deviceGenAndTransmitML(dev,am15,[]);
[pts,status]=sweepJVML(dev,G,0,1.60,0.02);
if ~status.completed, error('MD421 JV failed: %s',status.stopReason); end
sim=struct('pts',pts,'status',status,'metrics',pvMetrics(pts,spectrumPower(am15)));
end

function [bestRs,bestRsh,bestM,bestLoss]=fitTerminal(pts,Pin,target)
RsGrid=unique([0 logspace(-1,3,31)]);
RshGrid=logspace(1,7,31);
bestLoss=Inf; bestRs=NaN; bestRsh=NaN; bestM=struct();
for Rs=RsGrid
 for Rsh=RshGrid
  if Rs>=Rsh,continue;end
  m=pvMetrics(applyTerminalResistance(pts,Rs,Rsh),Pin);
  if ~m.valid,continue;end
  % Jsc/Voc are already bulk/optical gates. Terminal fitting emphasizes FF
  % while penalizing damage to those previously calibrated quantities.
  loss=((m.Jsc-target.Jsc)/(0.10*target.Jsc))^2 + ...
       ((m.Voc-target.Voc)/0.05)^2 + ...
       ((m.FF-target.FF)/0.05)^2;
  if loss<bestLoss
   bestLoss=loss; bestRs=Rs; bestRsh=Rsh; bestM=m;
  end
 end
end
end

function L=layer(t,N,Eg,chi,eps_,Nd,Na,muN,muP,tau,Nc,Nv,useTable,opt,centers)
L=struct('t',t,'N',N,'Eg',Eg,'chi',chi,'eps',eps_,'Nd',Nd,'Na',Na, ...
 'muN',muN,'muP',muP,'tau',tau,'Nc',Nc,'Nv',Nv,'NcNv',sqrt(Nc*Nv), ...
 'Brad',0,'alphaMax',0,'edgeWidth',0.05,'useOpticalTable',useTable, ...
 'opticalTable',opt,'srhCenters',centers);
end
