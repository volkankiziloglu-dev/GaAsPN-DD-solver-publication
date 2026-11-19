function T = run_paper_referee_checks(part)
%RUN_PAPER_REFEREE_CHECKS Two checks raised by the simulated referee review (docs/HAKEM_SIMULASYONU_2026-11-04.md).
% part 'A': n-/p-type absorber lifetime ratio. RHJ has an n-type absorber with lifetime r*tauP; HJ-0, F-HJ and RBSF-GaP keep tauP.
%           tauP in {600, 100, 10, 1} ns, r in {1, 3, 10, 30, 100}; AM1.5G and 6500 K / 1000 lx. Other parameters nominal.
% part 'B': rear-contact boundary of the HJ-0 control (and RBSF-GaP for reference): Sback in {1e5, 1e3, 0} cm/s.
% Output: results/paper_referee_checks/referee_checks_{A,B}.csv (checkpointed after every row).
if nargin < 1, part = 'A'; end
assertCanonicalPath(); amMode = strcmp(getenv('AM_UNIFIED'),'1'); icdRun = 2 - amMode; amSuffix = ''; if amMode, amSuffix = '_am'; end   % AM_UNIFIED=1: rerun AM1.5G with all finite points (unified extractor)
warning('off','MATLAB:nearlySingularMatrix'); warning('off','MATLAB:singularMatrix');
root = fileparts(fileparts(mfilename('fullpath'))); outDir = fullfile(root,'results',['paper_referee_checks_unified' amSuffix]); if ~isfolder(outDir), mkdir(outDir); end
arc = load(fullfile(root,'results','paper_stageA','stageA_comparison.mat'),'d1nm','d2nm');
arcFn = @(lam) tmmReflectance(lam,[1.0,1.38,2.35,gapRefractiveIndexReal(lam)],[arc.d1nm,arc.d2nm]);
Bmat = paperRadiativeCoefficients(300);
conds = {'AM1.5G'; '6500K_1000lx'};
specs = {loadRealAM15G(); buildCIEIndoorSpectrum('6500K',1000)};
Vgrids = {(0:0.02:1.70)'; (0:0.02:1.56)'};
Pins = cellfun(@spectrumPower,specs);
common = {'BradAbs',Bmat.GaAsPN_cm3ps,'BradGaP',Bmat.GaP_cm3ps,'BradGaNP',Bmat.GaPN_Ege18_cm3ps, ...
  'CnAugerAbs',3e-30,'CpAugerAbs',3e-30,'CnAugerGaP',3e-30,'CpAugerGaP',3e-30,'CnAugerGaNP',3e-30,'CpAugerGaNP',3e-30};
rows = {};
switch upper(part)
  case 'A'
    taus = [600 100 10 1]*1e-9; ratios = [1 3 10 30 100];
    jobs = {};
    for it = 1:numel(taus)
      for code = ["HJ-0","F-HJ","RBSF-GaP"]
        jobs(end+1,:) = {char(code), taus(it), 1, 1e5}; %#ok<AGROW>
      end
      for ir = 1:numel(ratios)
        jobs(end+1,:) = {'RHJ', taus(it), ratios(ir), 1e5}; %#ok<AGROW>
      end
    end
  case 'B'
    jobs = {};
    for S = [1e5 1e3 0]
      for code = ["HJ-0","RBSF-GaP"], jobs(end+1,:) = {char(code), 0.6e-6, 1, S}; end %#ok<AGROW>
    end
  case 'C'
    jobs = {};
    for S = [1e5 1e3 0], jobs(end+1,:) = {'F-HJ', 0.6e-6, 1, S}; end %#ok<AGROW>
  otherwise
    error('part must be A, B or C');
end
ck = fullfile(outDir,sprintf('referee_checks_%s.csv',upper(part)));
for j = 1:size(jobs,1)
  code = jobs{j,1}; tauP = jobs{j,2}; r = jobs{j,3}; Sb = jobs{j,4};
  tauAbs = tauP; if strcmp(code,'RHJ'), tauAbs = tauP*r; end
  for icd = icdRun   % indoor only (AM1.5G unchanged); unified extractor
    t0 = tic;
    try
      dev = buildPaperArchitecture(code, 'AbsorberTau_s',tauAbs,'Sback',Sb, common{:});
      G = deviceGenAndTransmitML(dev,specs{icd},arcFn);
      if icd == 1
        [~,~,m] = runAcceptedJVSweep(dev,G,Vgrids{icd},Pins(icd),true); cls = 'AM15G';
        Jsc = m.Jsc*1e3; Voc = m.Voc*1e3; FF = m.FF*100; PCE = m.PCE;
      else
        [~,gate] = indoorSweepLocal(dev,G,Vgrids{icd},Pins(icd)); cls = char(gate.curveClass);
        Jsc = gate.Jsc*1e6; Voc = gate.Voc*1e3; FF = gate.FF*100; PCE = gate.PCE;
      end
    catch ME
      Jsc = NaN; Voc = NaN; FF = NaN; PCE = NaN; cls = ['ERR ' strrep(ME.message,',',';')];
    end
    rows(end+1,:) = {code, tauP, r, tauAbs, Sb, conds{icd}, Jsc, Voc, FF, PCE, cls}; %#ok<AGROW>
    T = cell2table(rows,'VariableNames',{'Architecture','tauP_s','ratio','tauAbs_s','Sback','Condition','Jsc','Voc_mV','FF_pct','PCE_pct','Class'});
    writetable(T,ck);
    fprintf('%s tauAbs=%.1e S=%g %s: PCE=%.3f class=%s (%.0f s)\n',code,tauAbs,Sb,conds{icd},PCE,cls,toc(t0));
  end
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
