function T = run_paper_absorber_doping_ranking()
%RUN_PAPER_ABSORBER_DOPING_RANKING Does the architecture ranking survive a change of the GaAsPN absorber doping?
% The four main-text architectures are rebuilt with absorber/base doping N in {3e16, 1e17 (frozen default), 3e17} cm-3 (p-type for HJ-0, F-HJ and
% RBSF-GaP, n-type for RHJ; the HJ-0 n-emitter and the p+ layers are unchanged) and solved under AM1.5G and CIE 6500 K / 1000 lx with the production
% sweeps and acceptance gates (run_paper_stageA_comparison.m / run_paper_indoor_production.m). Full recombination (SRH + radiative + Auger).
% Two readings: (i) common doping change of all architectures; (ii) only the RHJ doping changed while the others stay at 1e17.
% Output: results/paper_absorber_doping/doping_ranking.csv.
assertCanonicalPath(); amMode = strcmp(getenv('AM_UNIFIED'),'1'); icdRun = 2 - amMode; amSuffix = ''; if amMode, amSuffix = '_am'; end   % AM_UNIFIED=1: rerun AM1.5G with all finite points (unified extractor)
warning('off','MATLAB:nearlySingularMatrix'); warning('off','MATLAB:singularMatrix');
root = fileparts(fileparts(mfilename('fullpath'))); outDir = fullfile(root,'results',['paper_absorber_doping_unified' amSuffix]); if ~isfolder(outDir), mkdir(outDir); end
arc = load(fullfile(root,'results','paper_stageA','stageA_comparison.mat'),'d1nm','d2nm');
arcFn = @(lam) tmmReflectance(lam,[1.0,1.38,2.35,gapRefractiveIndexReal(lam)],[arc.d1nm,arc.d2nm]);
Bmat = paperRadiativeCoefficients(300);
codes = ["HJ-0","F-HJ","RBSF-GaP","RHJ"]; Ns = [3e16 1e17 3e17];
conds = {'AM1.5G', 'AM15G'; '6500K_1000lx', '6500K'};
rows = {};
for icd = icdRun   % indoor only (AM1.5G unchanged); unified extractor
  if icd == 1, spec = loadRealAM15G(); Vgrid = (0:0.02:1.70)'; else, spec = buildCIEIndoorSpectrum('6500K',1000); Vgrid = (0:0.02:1.56)'; end
  Pin = spectrumPower(spec);
  for iN = 1:numel(Ns)
    for ia = 1:numel(codes)
      [dev,~] = buildPaperArchitecture(codes(ia), AbsorberDoping=Ns(iN), BradAbs=Bmat.GaAsPN_cm3ps, BradGaP=Bmat.GaP_cm3ps, BradGaNP=Bmat.GaPN_Ege18_cm3ps, ...
        CnAugerAbs=3e-30,CpAugerAbs=3e-30,CnAugerGaP=3e-30,CpAugerGaP=3e-30,CnAugerGaNP=3e-30,CpAugerGaNP=3e-30);
      G = deviceGenAndTransmitML(dev,spec,arcFn);
      try
        if icd == 1
          [~,~,m] = runAcceptedJVSweep(dev,G,Vgrid,Pin,true); cls = 'AM15G';
          Jsc = m.Jsc*1e3; Voc = m.Voc*1e3; FF = m.FF*100; PCE = m.PCE;
        else
          [~,gate] = indoorSweepLocal(dev,G,Vgrid,Pin); cls = char(gate.curveClass);
          Jsc = gate.Jsc*1e6; Voc = gate.Voc*1e3; FF = gate.FF*100; PCE = gate.PCE;
        end
      catch ME
        Jsc = NaN; Voc = NaN; FF = NaN; PCE = NaN; cls = ['ERR ' ME.message]; %#ok<NASGU>
      end
      fprintf('%-13s N=%.0e %-9s: Jsc %.4g, Voc %.1f mV, FF %.2f, PCE %.3f (%s)\n', conds{icd,1}, Ns(iN), codes(ia), Jsc, Voc, FF, PCE, cls);
      rows(end+1,:) = {conds{icd,1}, Ns(iN), char(codes(ia)), Jsc, Voc, FF, PCE, cls}; %#ok<AGROW>
    end
  end
end
T = cell2table(rows,'VariableNames',{'Condition','N_cm3','Architecture','Jsc','Voc_mV','FF_pct','PCE_pct','Class'});
writetable(T,fullfile(outDir,'doping_ranking.csv'));
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
