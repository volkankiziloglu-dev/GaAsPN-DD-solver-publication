function T = run_paper_referee_checks2(part)
%RUN_PAPER_REFEREE_CHECKS2 Further robustness checks (absorption edge and layer thickness).
% part 'E': GaAsPN optical-edge shift {-100,-50,0,+50,+100} meV (rigid shift of the optical table; Eg of the electrical model unchanged), all four architectures.
% part 'T': absorber thickness scale {1,2} x rear-contact velocity {1e5,1e3,0} cm/s, all four architectures (HJ-0 base 1350 -> 2700 nm; other absorbers 1500 -> 3000 nm etc.).
% AM1.5G and CIE 6500 K / 1000 lx, production sweeps and acceptance gates, nominal parameters otherwise.
% Output: results/paper_referee_checks/referee_checks_<part>.csv.
assertCanonicalPath(); amMode = strcmp(getenv('AM_UNIFIED'),'1'); icdRun = 2 - amMode; amSuffix = ''; if amMode, amSuffix = '_am'; end   % AM_UNIFIED=1: rerun AM1.5G with all finite points (unified extractor) warning('off','all');
root = fileparts(fileparts(mfilename('fullpath'))); outDir = fullfile(root,'results',['paper_referee_checks_unified' amSuffix]); if ~isfolder(outDir), mkdir(outDir); end
arc = load(fullfile(root,'results','paper_stageA','stageA_comparison.mat'),'d1nm','d2nm');
arcFn = @(lam) tmmReflectance(lam,[1.0,1.38,2.35,gapRefractiveIndexReal(lam)],[arc.d1nm,arc.d2nm]);
Bmat = paperRadiativeCoefficients(300);
common = {'BradAbs',Bmat.GaAsPN_cm3ps,'BradGaP',Bmat.GaP_cm3ps,'BradGaNP',Bmat.GaPN_Ege18_cm3ps, ...
  'CnAugerAbs',3e-30,'CpAugerAbs',3e-30,'CnAugerGaP',3e-30,'CpAugerGaP',3e-30,'CnAugerGaNP',3e-30,'CpAugerGaNP',3e-30};
conds = {'AM1.5G'; '6500K_1000lx'};
specs = {loadRealAM15G(); buildCIEIndoorSpectrum('6500K',1000)};
Vgrids = {(0:0.02:1.70)'; (0:0.02:1.56)'};
Pins = cellfun(@spectrumPower,specs);
codes = {'HJ-0','F-HJ','RBSF-GaP','RHJ'};
jobs = {};   % {code, shift_eV, thicknessScale, Sback}
switch upper(part)
  case 'E'
    for sh = [-0.10 -0.05 0 0.05 0.10], for ic = 1:4, jobs(end+1,:) = {codes{ic}, sh, 1, 1e5}; end, end %#ok<AGROW>
  case 'T'
    for ts = [1 2], for S = [1e5 1e3 0], for ic = 1:4, jobs(end+1,:) = {codes{ic}, 0, ts, S}; end, end, end %#ok<AGROW>
  otherwise, error('part must be E or T');
end
rows = {}; ck = fullfile(outDir,sprintf('referee_checks_%s.csv',upper(part)));
for j = 1:size(jobs,1)
  code = jobs{j,1}; sh = jobs{j,2}; ts = jobs{j,3}; Sb = jobs{j,4};
  for icd = icdRun   % indoor only (AM1.5G unchanged); unified extractor
    t0 = tic;
    try
      dev = buildPaperArchitecture(code,'AbsorberOpticalShift_eV',sh,'AbsorberThicknessScale',ts,'Sback',Sb,common{:});
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
    rows(end+1,:) = {code, sh, ts, Sb, conds{icd}, Jsc, Voc, FF, PCE, cls}; %#ok<AGROW>
    T = cell2table(rows,'VariableNames',{'Architecture','shift_eV','thickScale','Sback','Condition','Jsc','Voc_mV','FF_pct','PCE_pct','Class'});
    writetable(T,ck);
    fprintf('%s shift=%+.2f ts=%g S=%g %s: PCE=%.3f %s (%.0f s)\n',code,sh,ts,Sb,conds{icd},PCE,cls,toc(t0));
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
