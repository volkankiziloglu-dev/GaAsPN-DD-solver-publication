function T = run_paper_pin(part,slice)
%RUN_PAPER_PIN n+-GaP / i-GaAsPN / p+-GaP (PIN) architecture, the structure of the grown GaAsPN cells (simulated referees R2-4, R5-M2; docs/HAKEM_SIMULASYONU_2026-11-06.md).
% Same optics (ARC, TMM generation), recombination (SRH + B_rad + Auger 3e-30), contacts and unified extractor (all finite points, J_sc at 0 V,
% linear V_oc, parabolic P_max) as the lifetime table of the paper (run_paper_referee_checks part A, merged results).
% part 'N': PIN at the common 1.5 um GaAsPN thickness; i-layer background +1e15 (n), -1e15 (p), +1e16 (n) cm-3; tau in {600,100,10,1} ns
%           (the i-layer has one lifetime for both carriers, tau_n = tau_p = tau; ratio r is not defined for it).
% part 'M': mu_n = 50 cm2/Vs; PIN (p-type 1e15) and the four architectures at tau_p in {600,10,1} ns, 1.5 um.
% part 'U' / 'V': HJ-U (p+-GaAsPN rear layer) control, mu_n = 500 (tau 600/100/10/1 ns) / 50 (600/10/1 ns).
% part 'T': thin absorbers (thickness scale 1/3, 0.5 um GaAsPN): PIN (+1e15) and HJ-0, F-HJ, RBSF-GaP, RHJ (r = 1 and 10) at tau_p in {600,10,1} ns.
% Output: results/paper_pin/pin_<part>_<k>of<n>.csv (checkpointed after every row).
if nargin < 1, part = 'N'; end
assertCanonicalPath();
warning('off','MATLAB:nearlySingularMatrix'); warning('off','MATLAB:singularMatrix');
root = fileparts(fileparts(mfilename('fullpath'))); outDir = fullfile(root,'results','paper_pin'); if ~isfolder(outDir), mkdir(outDir); end
arc = load(fullfile(root,'results','paper_stageA','stageA_comparison.mat'),'d1nm','d2nm');
arcFn = @(lam) tmmReflectance(lam,[1.0,1.38,2.35,gapRefractiveIndexReal(lam)],[arc.d1nm,arc.d2nm]);
Bmat = paperRadiativeCoefficients(300);
conds = {'AM1.5G'; '6500K_1000lx'};
specs = {loadRealAM15G(); buildCIEIndoorSpectrum('6500K',1000)};
Vgrids = {(0:0.02:1.70)'; (0:0.02:1.56)'};
Pins = cellfun(@spectrumPower,specs);
common = {'BradAbs',Bmat.GaAsPN_cm3ps,'BradGaP',Bmat.GaP_cm3ps,'BradGaNP',Bmat.GaPN_Ege18_cm3ps, ...
  'CnAugerAbs',3e-30,'CpAugerAbs',3e-30,'CnAugerGaP',3e-30,'CpAugerGaP',3e-30,'CnAugerGaNP',3e-30,'CpAugerGaNP',3e-30};
% job = {code, tauP, ratio, background, thicknessScale}
jobs = {};
switch upper(part)
  case 'N'
    for tau = [600 100 10 1]*1e-9
      for bg = [1e15 -1e15 1e16], jobs(end+1,:) = {'PIN', tau, 1, bg, 1}; end %#ok<AGROW>
    end
  case 'T'
    for tau = [600 10 1]*1e-9
      jobs(end+1,:) = {'PIN', tau, 1, 1e15, 1/3}; %#ok<AGROW>
      for code = ["HJ-0","F-HJ","RBSF-GaP"], jobs(end+1,:) = {char(code), tau, 1, NaN, 1/3}; end %#ok<AGROW>
      for r = [1 10], jobs(end+1,:) = {'RHJ', tau, r, NaN, 1/3}; end %#ok<AGROW>
    end
  case 'M'
    % low electron mobility (mu_n = 50 cm2/Vs instead of 500): PIN (p-type 1e15 background) and the four architectures at tau_p in {600,10,1} ns, 1.5 um
    for tau = [600 10 1]*1e-9
      jobs(end+1,:) = {'PIN', tau, 1, -1e15, 1}; %#ok<AGROW>
      for code = ["HJ-0","F-HJ","RBSF-GaP"], jobs(end+1,:) = {char(code), tau, 1, NaN, 1}; end %#ok<AGROW>
      for r = [1 10], jobs(end+1,:) = {'RHJ', tau, r, NaN, 1}; end %#ok<AGROW>
    end
  case 'U'
    % homojunction with p+-GaAsPN rear layer (HJ-U) as the passivated-rear homojunction control (referees R2-4 / R5-M1), tau_p in {600,100,10,1} ns
    for tau = [600 100 10 1]*1e-9, jobs(end+1,:) = {'HJ-U', tau, 1, NaN, 1}; end %#ok<AGROW>
  case 'V'
    % HJ-U with mu_n = 50 cm2/Vs at tau_p in {600,10,1} ns (companion of part M)
    for tau = [600 10 1]*1e-9, jobs(end+1,:) = {'HJ-U', tau, 1, NaN, 1}; end %#ok<AGROW>
  otherwise
    error('part must be N, T, M, U or V');
end
muN = 500; if any(upper(part) == ['M' 'V']), muN = 50; end
% slice = [k n]: run jobs k, k+n, k+2n, ... (parallel workers); output pin_<part>_<k>of<n>.csv
if nargin < 2 || isempty(slice), slice = [1 1]; end
jobs = jobs(slice(1):slice(2):end,:);
rows = {}; ck = fullfile(outDir,sprintf('pin_%s_%dof%d.csv',upper(part),slice(1),slice(2)));
for j = 1:size(jobs,1)
  [code,tauP,r,bg,tS] = jobs{j,:};
  tauAbs = tauP; if strcmp(code,'RHJ'), tauAbs = tauP*r; end
  extra = {}; if strcmp(code,'PIN'), extra = {'PinBackground',bg}; end
  for icd = 1:2
    t0 = tic;
    try
      dev = buildPaperArchitecture(code,'AbsorberTau_s',tauAbs,'AbsorberThicknessScale',tS,'AbsorberMuN',muN,extra{:},common{:});
      G = deviceGenAndTransmitML(dev,specs{icd},arcFn);
      if icd == 1
        [~,~,m] = runAcceptedJVSweep(dev,G,Vgrids{icd},Pins(icd),true); cls = 'AM15G';
        Jsc = m.Jsc*1e3; Voc = m.Voc*1e3; FF = m.FF*100; PCE = m.PCE;
      else
        gate = indoorSweepLocal(dev,G,Vgrids{icd},Pins(icd)); cls = char(gate.curveClass);
        Jsc = gate.Jsc*1e6; Voc = gate.Voc*1e3; FF = gate.FF*100; PCE = gate.PCE;
      end
    catch ME
      Jsc = NaN; Voc = NaN; FF = NaN; PCE = NaN; cls = ['ERR ' strrep(ME.message,',',';')];
    end
    rows(end+1,:) = {code, tauP, r, tauAbs, bg, tS, muN, conds{icd}, Jsc, Voc, FF, PCE, cls}; %#ok<AGROW>
    T = cell2table(rows,'VariableNames',{'Architecture','tauP_s','ratio','tauAbs_s','PinBackground','ThicknessScale','MuN','Condition','Jsc','Voc_mV','FF_pct','PCE_pct','Class'});
    writetable(T,ck);
    fprintf('%s tau=%.1e bg=%g tS=%.3f %s: PCE=%.3f class=%s (%.0f s)\n',code,tauAbs,bg,tS,conds{icd},PCE,cls,toc(t0));
  end
end
end

function gate = indoorSweepLocal(dev,G,Vgrid,Pin)
% Identical to run_paper_referee_checks>indoorSweepLocal (unified extractor).
N = dev.N; idxStart = find(abs(Vgrid-0.02)<1e-12,1);
r002 = gummelSolveQFPAnalytic(dev,0.02,G,struct('psiInit',dev.psiEq,'phinInit',zeros(N,1),'phipInit',zeros(N,1)));
prev = r002;
for v = [0.01,0.005,0], r0 = gummelSolveQFPAnalytic(dev,v,G,struct('psiInit',prev.psi,'phinInit',prev.phin,'phipInit',prev.phip)); prev = r0; end
states = repmat(r002,numel(Vgrid),1); states(1) = r0; states(idxStart) = r002; prev = r002;
for k = idxStart+1:numel(Vgrid)
  rk = gummelSolveQFPAnalytic(dev,Vgrid(k),G,struct('psiInit',prev.psi,'phinInit',prev.phin,'phipInit',prev.phip));
  states(k) = rk; if rk.isFinite, prev = rk; end
end
gate = indoorAcceptanceGate(dev,G,Vgrid,states,Pin);
fU = gate.isFinite(:); mU = pvMetrics([gate.Vs(fU), gate.Jterm(fU)],Pin);
gate.Jsc = mU.Jsc; gate.Voc = mU.Voc; gate.FF = mU.FF; gate.PCE = mU.PCE;
end
