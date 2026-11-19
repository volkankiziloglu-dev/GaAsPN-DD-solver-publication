function T = run_paper_sfront(slice)
%RUN_PAPER_SFRONT Finite front-contact recombination velocity (simulated referees R5-M4, R1; docs/HAKEM_SIMULASYONU_2026-11-06.md).
% The production model uses an ideal ohmic (Dirichlet) front contact, i.e. S_front = infinity for the minority carrier (holes in the
% n+-GaP window). With the optional field dev.Sfront the minority-carrier front row becomes the Robin law -J/q = S_front*(p - p0)
% (assembleResidualQFP.m and the analytic Jacobians, verified by cases/check_sfront_impl.m; without the field the production result
% is unchanged). Production sweeps and the unified extractor (all finite points, J_sc at 0 V, linear V_oc, parabolic P_max).
% Jobs: HJ-0, F-HJ, RBSF-GaP, RHJ (r = 1) x S_front in {Inf, 1e6, 1e4, 1e2} cm/s x tau_p in {600, 10, 1} ns x {AM1.5G, 6500 K / 1000 lx}.
% Output: results/paper_sfront/sfront_<k>of<n>.csv.
if nargin < 1 || isempty(slice), slice = [1 1]; end
assertCanonicalPath(); warning('off','all');
root = fileparts(fileparts(mfilename('fullpath'))); outDir = fullfile(root,'results','paper_sfront'); if ~isfolder(outDir), mkdir(outDir); end
arc = load(fullfile(root,'results','paper_stageA','stageA_comparison.mat'),'d1nm','d2nm');
arcFn = @(lam) tmmReflectance(lam,[1.0,1.38,2.35,gapRefractiveIndexReal(lam)],[arc.d1nm,arc.d2nm]);
Bmat = paperRadiativeCoefficients(300);
conds = {'AM1.5G'; '6500K_1000lx'};
specs = {loadRealAM15G(); buildCIEIndoorSpectrum('6500K',1000)};
Vgrids = {(0:0.02:1.70)'; (0:0.02:1.56)'};
Pins = cellfun(@spectrumPower,specs);
common = {'BradAbs',Bmat.GaAsPN_cm3ps,'BradGaP',Bmat.GaP_cm3ps,'BradGaNP',Bmat.GaPN_Ege18_cm3ps, ...
  'CnAugerAbs',3e-30,'CpAugerAbs',3e-30,'CnAugerGaP',3e-30,'CpAugerGaP',3e-30,'CnAugerGaNP',3e-30,'CpAugerGaNP',3e-30};
jobs = {};
for tau = [600 10 1]*1e-9
  for Sf = [Inf 1e6 1e4 1e2]
    for code = ["HJ-0","F-HJ","RBSF-GaP","RHJ"], jobs(end+1,:) = {char(code), tau, Sf}; end %#ok<AGROW>
  end
end
jobs = jobs(slice(1):slice(2):end,:);
rows = {}; ck = fullfile(outDir,sprintf('sfront_%dof%d.csv',slice(1),slice(2)));
for j = 1:size(jobs,1)
  [code,tauP,Sf] = jobs{j,:};
  for icd = 1:2
    t0 = tic;
    try
      dev = buildPaperArchitecture(code,'AbsorberTau_s',tauP,common{:});
      if isfinite(Sf), dev.Sfront = Sf; end
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
    rows(end+1,:) = {code, tauP, Sf, conds{icd}, Jsc, Voc, FF, PCE, cls}; %#ok<AGROW>
    T = cell2table(rows,'VariableNames',{'Architecture','tauP_s','Sfront','Condition','Jsc','Voc_mV','FF_pct','PCE_pct','Class'});
    writetable(T,ck);
    fprintf('%s tau=%.1e Sf=%g %s: PCE=%.4f %s (%.0f s)\n',code,tauP,Sf,conds{icd},PCE,cls,toc(t0));
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
