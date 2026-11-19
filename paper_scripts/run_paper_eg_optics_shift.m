function T = run_paper_eg_optics_shift(variant)
%RUN_PAPER_EG_OPTICS_SHIFT Referee check: absorption edge shifted TOGETHER with the electrical bandgap (physically consistent version of the optical-edge test).
% The GaAsPN optical table is shifted rigidly by dE in {-100,-50,0,+50,+100} meV, the electrical bandgap of the absorber is shifted by the same amount
% (so that n_i changes) and the radiative coefficient is recomputed with the van Roosbroeck-Shockley integral from the shifted table and the shifted gap.
% variant 'CB': the shift is taken in the conduction band (chi_abs = 4.0 - dE, valence band fixed; nitrogen acts on the conduction band);
% variant 'VB': electron affinity fixed (the valence band moves). All four architectures, AM1.5G and CIE 6500 K / 1000 lx, production sweeps and gates, nominal parameters otherwise.
% Output: results/paper_referee_checks/eg_optics_shift_<variant>.csv.
assertCanonicalPath(); amMode = strcmp(getenv('AM_UNIFIED'),'1'); icdRun = 2 - amMode; amSuffix = ''; if amMode, amSuffix = '_am'; end   % AM_UNIFIED=1: rerun AM1.5G with all finite points (unified extractor) warning('off','all');
root = fileparts(fileparts(mfilename('fullpath'))); outDir = fullfile(root,'results',['paper_referee_checks_unified' amSuffix]); if ~isfolder(outDir), mkdir(outDir); end
arc = load(fullfile(root,'results','paper_stageA','stageA_comparison.mat'),'d1nm','d2nm');
arcFn = @(lam) tmmReflectance(lam,[1.0,1.38,2.35,gapRefractiveIndexReal(lam)],[arc.d1nm,arc.d2nm]);
Bmat = paperRadiativeCoefficients(300);
absBand = gaasnpBandgap(0.03);
optAbs = loadGaAsNPOptical(3.0);
conds = {'AM1.5G'; '6500K_1000lx'};
specs = {loadRealAM15G(); buildCIEIndoorSpectrum('6500K',1000)};
Vgrids = {(0:0.02:1.70)'; (0:0.02:1.56)'};
Pins = cellfun(@spectrumPower,specs);
codes = {'HJ-0','F-HJ','RBSF-GaP','RHJ'};
inCB = strcmpi(variant,'CB');
rows = {}; ck = fullfile(outDir,sprintf('eg_optics_shift_%s.csv',upper(variant)));
for sh = [-0.10 -0.05 0 0.05 0.10]
  o = optAbs; o.E = o.E + sh;
  Bsh = computeBradVRS(o,absBand.Eg + sh,1e19,1e19,300);
  common = {'BradAbs',Bsh,'BradGaP',Bmat.GaP_cm3ps,'BradGaNP',Bmat.GaPN_Ege18_cm3ps, ...
    'CnAugerAbs',3e-30,'CpAugerAbs',3e-30,'CnAugerGaP',3e-30,'CpAugerGaP',3e-30,'CnAugerGaNP',3e-30,'CpAugerGaNP',3e-30};
  for ic = 1:4
    code = codes{ic};
    for icd = icdRun   % indoor only (AM1.5G unchanged); unified extractor
      t0 = tic;
      try
        dev = buildPaperArchitecture(code,'AbsorberOpticalShift_eV',sh,'AbsorberEgShift_eV',sh,'AbsorberEgShiftInCB',inCB,common{:});
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
      rows(end+1,:) = {code, sh, upper(variant), Bsh, conds{icd}, Jsc, Voc, FF, PCE, cls}; %#ok<AGROW>
      T = cell2table(rows,'VariableNames',{'Architecture','shift_eV','variant','BradAbs','Condition','Jsc','Voc_mV','FF_pct','PCE_pct','Class'});
      writetable(T,ck);
      fprintf('%s %s shift=%+.2f %s: PCE=%.3f %s (%.0f s)',variant,code,sh,conds{icd},PCE,cls,toc(t0)); disp(' ');
    end
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
