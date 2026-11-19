function T = run_paper_uncertainty_sampling(nSamples, nStart)
%RUN_PAPER_UNCERTAINTY_SAMPLING Paired parameter-uncertainty sampling of the four main architectures.
% Every sample draws one parameter set (Latin hypercube, fixed seed) and applies it to HJ-0, F-HJ, RBSF-GaP and RHJ, so architectures are compared
% pairwise at equal parameters. Samples are solved under AM1.5G and CIE 6500 K / 1000 lx with the production sweeps and acceptance gates; full
% recombination (SRH + radiative + Auger). Not a calibration: ranges are plausibility ranges around the frozen nominal values, not fits.
% Sampled (independent, ranges in parentheses):
%   tauAbs   GaAsPN SRH lifetime, log-uniform [0.1, 3] us (nominal 0.6 us)
%   Sint     GaP/GaAsPN interface velocity (front and rear), log-uniform [1e3, 1e5] cm/s (nominal 8910)
%   Nabs     absorber/base doping, log-uniform [3e16, 3e17] cm-3 (nominal 1e17)
%   chiGaP   GaP electron affinity, uniform [3.7, 4.3] eV (nominal 3.8; sets Delta E_C, Delta E_V of the GaP/GaAsPN heterointerface)
%   tauGaP   GaP minority lifetime, log-uniform [1e-11, 1e-8] s (nominal 1e-9)
%   muN      GaAsPN electron mobility, log-uniform [250, 1000] cm2/Vs (nominal 500)
% Output: results/paper_uncertainty/uncertainty_samples.csv (checkpointed after every sample).
if nargin < 1, nSamples = 60; end
if nargin < 2, nStart = 1; end
assertCanonicalPath(); amMode = strcmp(getenv('AM_UNIFIED'),'1'); icdRun = 2 - amMode; amSuffix = ''; if amMode, amSuffix = '_am'; end   % AM_UNIFIED=1: rerun AM1.5G with all finite points (unified extractor)
warning('off','MATLAB:nearlySingularMatrix'); warning('off','MATLAB:singularMatrix');
root = fileparts(fileparts(mfilename('fullpath'))); outDir = fullfile(root,'results',['paper_uncertainty_unified' amSuffix]); if ~isfolder(outDir), mkdir(outDir); end
arc = load(fullfile(root,'results','paper_stageA','stageA_comparison.mat'),'d1nm','d2nm');
arcFn = @(lam) tmmReflectance(lam,[1.0,1.38,2.35,gapRefractiveIndexReal(lam)],[arc.d1nm,arc.d2nm]);
Bmat = paperRadiativeCoefficients(300);
codes = ["HJ-0","F-HJ","RBSF-GaP","RHJ"];
conds = {'AM1.5G'; '6500K_1000lx'};
specs = {loadRealAM15G(); buildCIEIndoorSpectrum('6500K',1000)};
Vgrids = {(0:0.02:1.70)'; (0:0.02:1.56)'};
Pins = cellfun(@spectrumPower,specs);

% Latin hypercube in [0,1]^6, fixed seed (the design depends on nSamples: randperm(nSamples) changes with it, so a run with another nSamples is a different design)
rng(20261030,'twister'); d = 6; U = zeros(nSamples,d);
for j = 1:d, U(:,j) = (randperm(nSamples)' - rand(nSamples,1))/nSamples; end
lg = @(u,a,b) 10.^(log10(a) + u*(log10(b)-log10(a)));
P = table(lg(U(:,1),1e-7,3e-6), lg(U(:,2),1e3,1e5), lg(U(:,3),3e16,3e17), 3.7+0.6*U(:,4), lg(U(:,5),1e-11,1e-8), lg(U(:,6),250,1000), ...
  'VariableNames',{'tauAbs','Sint','Nabs','chiGaP','tauGaP','muN'});

nEnd = nSamples; if ~isempty(getenv('AM_END')), nEnd = str2double(getenv('AM_END')); end
ckName = 'uncertainty_samples.csv'; if amMode, ckName = sprintf('uncertainty_samples_%d_%d.csv',nStart,nEnd); end
ck = fullfile(outDir,ckName); rows = {};
if nStart > 1 && isfile(ck), old = readtable(ck); rows = table2cell(old); end
for is = nStart:nEnd
  t0 = tic;
  for icd = icdRun   % indoor only (AM1.5G unchanged); unified extractor
    spec = specs{icd}; Vgrid = Vgrids{icd}; Pin = Pins(icd);
    for ia = 1:numel(codes)
      try
        dev = buildPaperArchitecture(codes(ia), AbsorberDoping=P.Nabs(is), AbsorberTau_s=P.tauAbs(is), AbsorberMuN=P.muN(is), ...
          frontIRV=P.Sint(is), rearIRV=P.Sint(is), GaPChi_eV=P.chiGaP(is), GaPTau_s=P.tauGaP(is), ...
          BradAbs=Bmat.GaAsPN_cm3ps, BradGaP=Bmat.GaP_cm3ps, BradGaNP=Bmat.GaPN_Ege18_cm3ps, ...
          CnAugerAbs=3e-30,CpAugerAbs=3e-30,CnAugerGaP=3e-30,CpAugerGaP=3e-30,CnAugerGaNP=3e-30,CpAugerGaNP=3e-30);
        G = deviceGenAndTransmitML(dev,spec,arcFn);
        if icd == 1
          [~,~,m] = runAcceptedJVSweep(dev,G,Vgrid,Pin,true); cls = 'AM15G';
          Jsc = m.Jsc*1e3; Voc = m.Voc*1e3; FF = m.FF*100; PCE = m.PCE;
        else
          [~,gate] = indoorSweepLocal(dev,G,Vgrid,Pin); cls = char(gate.curveClass);
          Jsc = gate.Jsc*1e6; Voc = gate.Voc*1e3; FF = gate.FF*100; PCE = gate.PCE;
        end
      catch ME
        Jsc = NaN; Voc = NaN; FF = NaN; PCE = NaN; cls = ['ERR ' strrep(ME.message,',',';')];
      end
      rows(end+1,:) = {is, P.tauAbs(is), P.Sint(is), P.Nabs(is), P.chiGaP(is), P.tauGaP(is), P.muN(is), conds{icd}, char(codes(ia)), Jsc, Voc, FF, PCE, cls}; %#ok<AGROW>
    end
  end
  T = cell2table(rows,'VariableNames',{'Sample','tauAbs','Sint','Nabs','chiGaP','tauGaP','muN','Condition','Architecture','Jsc','Voc_mV','FF_pct','PCE_pct','Class'});
  writetable(T,ck);
  fprintf('sample %d/%d done in %.1f s\n', is, nSamples, toc(t0));
end
end

function [curve,gate] = indoorSweepLocal(dev,G,Vgrid,Pin)
N = dev.N; idxStart = find(abs(Vgrid-0.02)<1e-12,1);
r002 = gummelSolveQFPAnalytic(dev,0.02,G,struct('psiInit',dev.psiEq,'phinInit',zeros(N,1),'phipInit',zeros(N,1)));
prev = r002;
for v = [0.01,0.005,0], r0 = gummelSolveQFPAnalytic(dev,v,G,struct('psiInit',prev.psi,'phinInit',prev.phin,'phipInit',prev.phip)); prev = r0; end
states = repmat(r002,numel(Vgrid),1); states(1) = r0; states(idxStart) = r002; prev = r002;
for k = idxStart+1:numel(Vgrid)
  rk = gummelSolveQFPAnalytic(dev,Vgrid(k),G,struct('psiInit',prev.psi,'phinInit',prev.phin,'phipInit',prev.phip));
  states(k) = rk; if rk.isFinite, prev = rk; end
end
gate = indoorAcceptanceGate(dev,G,Vgrid,states,Pin); curve = gate.Jterm;
fU = gate.isFinite(:); mU = pvMetrics([gate.Vs(fU), gate.Jterm(fU)],Pin);   % unified extractor: all finite points, J_sc at 0 V, linear V_oc, parabolic MPP
gate.Jsc = mU.Jsc; gate.Voc = mU.Voc; gate.FF = mU.FF; gate.PCE = mU.PCE;
end
