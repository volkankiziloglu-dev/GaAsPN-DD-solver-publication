function T = run_paper_design_map(tauList_ns, tag)
%RUN_PAPER_DESIGN_MAP (tau_p, r, S_back) design map of RBSF-GaP versus RHJ (and the HJ-0 reference), all other parameters nominal.
% For every p-type lifetime tau_p in tauList_ns (ns), RBSF-GaP and HJ-0 use tau_p, and RHJ uses r*tau_p (n-type absorber) for a list of ratios r.
% S_back = 1e5 and 1e3 cm/s; AM1.5G and CIE 6500 K / 1000 lx; production sweeps and acceptance gates as in run_paper_uncertainty_sampling_ns.m.
% Output: results/paper_design_map/design_map_<tag>.csv (checkpointed after every tau_p).
if nargin < 1, tauList_ns = [1 3 10 30 100 300 600]; end
if nargin < 2, tag = 'all'; end
assertCanonicalPath(); amMode = strcmp(getenv('AM_UNIFIED'),'1'); icdRun = 2 - amMode; amSuffix = ''; if amMode, amSuffix = '_am'; end   % AM_UNIFIED=1: rerun AM1.5G with all finite points (unified extractor)
warning('off','MATLAB:nearlySingularMatrix'); warning('off','MATLAB:singularMatrix');
root = fileparts(fileparts(mfilename('fullpath'))); outDir = fullfile(root,'results',['paper_design_map_unified' amSuffix]); if ~isfolder(outDir), mkdir(outDir); end
arc = load(fullfile(root,'results','paper_stageA','stageA_comparison.mat'),'d1nm','d2nm');
arcFn = @(lam) tmmReflectance(lam,[1.0,1.38,2.35,gapRefractiveIndexReal(lam)],[arc.d1nm,arc.d2nm]);
Bmat = paperRadiativeCoefficients(300);
conds = {'AM1.5G'; '6500K_1000lx'};
specs = {loadRealAM15G(); buildCIEIndoorSpectrum('6500K',1000)};
Vgrids = {(0:0.02:1.70)'; (0:0.02:1.56)'};
Pins = cellfun(@spectrumPower,specs);
rList = [1 1.5 2 3 5 10 20 30 60];
SbList = [1e5 1e3];
ck = fullfile(outDir,sprintf('design_map_%s.csv',tag)); rows = {};
for tp = tauList_ns(:)'
  t0 = tic;
  for Sb = SbList
    for icd = icdRun   % indoor only (AM1.5G unchanged); unified extractor
      spec = specs{icd}; Vgrid = Vgrids{icd}; Pin = Pins(icd);
      jobs = {'HJ-0',1; 'RBSF-GaP',1};
      for r = rList, jobs(end+1,:) = {'RHJ',r}; end %#ok<AGROW>
      for ij = 1:size(jobs,1)
        code = jobs{ij,1}; r = jobs{ij,2}; tauAbs = tp*1e-9*r;
        try
          dev = buildPaperArchitecture(code,'AbsorberTau_s',tauAbs,'Sback',Sb, ...
            'BradAbs',Bmat.GaAsPN_cm3ps,'BradGaP',Bmat.GaP_cm3ps,'BradGaNP',Bmat.GaPN_Ege18_cm3ps, ...
            'CnAugerAbs',3e-30,'CpAugerAbs',3e-30,'CnAugerGaP',3e-30,'CpAugerGaP',3e-30,'CnAugerGaNP',3e-30,'CpAugerGaNP',3e-30);
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
        rows(end+1,:) = {tp, r, Sb, conds{icd}, code, tauAbs, Jsc, Voc, FF, PCE, cls}; %#ok<AGROW>
      end
    end
  end
  T = cell2table(rows,'VariableNames',{'tauP_ns','ratio','Sback','Condition','Architecture','tauAbs_s','Jsc','Voc_mV','FF_pct','PCE_pct','Class'});
  writetable(T,ck);
  fprintf('tau_p = %g ns done in %.1f s\n', tp, toc(t0));
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
