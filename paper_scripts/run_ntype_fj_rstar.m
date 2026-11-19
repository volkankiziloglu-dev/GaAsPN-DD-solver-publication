function T = run_ntype_fj_rstar(part, nParts)
%RUN_NTYPE_FJ_RSTAR Crossover ratio r* of the n-type front-junction cell (RBSF-GaP-n: p+-GaP / n-GaAsPN / n+-GaP) against
% RBSF-GaP, on the grid of the design map (run_paper_design_map.m: tau_p = 1,3,10,30,100,300,600 ns, r = 1,1.5,2,3,5,10,20,30,60,
% S_back = 1e5 cm/s, AM1.5G and CIE 6500 K / 1000 lx, unified extractor). RBSF-GaP-n depends only on its own absorber lifetime
% tau_n = r*tau_p, so it is solved once per distinct tau_n; this run solves the share 'part' of 'nParts' of these lifetimes.
% Sweep: mirrored (V -> -V, J -> -J; p-type front contact), checked point by point with the native continuity residuals
% (see run_ntype_front_junction_check.m). r* is evaluated by analyze_ntype_fj_rstar.m.
% Output: results/paper_ntype_front_junction/ntype_fj_taun_part<k>.csv
if nargin < 1, part = 1; nParts = 1; end
assertCanonicalPath(); warning('off','all');
root = fileparts(fileparts(mfilename('fullpath'))); outDir = fullfile(root,'results','paper_ntype_front_junction'); if ~isfolder(outDir), mkdir(outDir); end
arc = load(fullfile(root,'results','paper_stageA','stageA_comparison.mat'),'d1nm','d2nm');
arcFn = @(lam) tmmReflectance(lam,[1.0,1.38,2.35,gapRefractiveIndexReal(lam)],[arc.d1nm,arc.d2nm]);
Bmat = paperRadiativeCoefficients(300);
conds = {'AM1.5G'; '6500K_1000lx'}; specs = {loadRealAM15G(); buildCIEIndoorSpectrum('6500K',1000)};
Vgrids = {(0:0.02:1.70)'; (0:0.02:1.56)'}; Pins = cellfun(@spectrumPower,specs);
tauP = [1 3 10 30 100 300 600]; rList = [1 1.5 2 3 5 10 20 30 60];
tauN = unique(round(tauP(:)*rList(:)',6)); tauN = tauN(part:nParts:end);
rows = {};
for tn = tauN(:)'
  dev = buildPaperArchitecture("RBSF-GaP-n",'AbsorberTau_s',tn*1e-9,'Sback',1e5, ...
    'BradAbs',Bmat.GaAsPN_cm3ps,'BradGaP',Bmat.GaP_cm3ps,'BradGaNP',Bmat.GaPN_Ege18_cm3ps, ...
    'CnAugerAbs',3e-30,'CpAugerAbs',3e-30,'CnAugerGaP',3e-30,'CpAugerGaP',3e-30,'CnAugerGaNP',3e-30,'CpAugerGaNP',3e-30);
  for icd = 1:2
    G = deviceGenAndTransmitML(dev,specs{icd},arcFn);
    [m,nOk,nRes,nV] = mirroredSweep(dev,G,Vgrids{icd},Pins(icd));
    sc = 1e3*(icd==1) + 1e6*(icd==2);
    rows(end+1,:) = {tn, conds{icd}, m.Jsc*sc, m.Voc*1e3, m.FF*100, m.PCE, nOk, nRes, nV}; %#ok<AGROW>
    fprintf('tau_n %8.1f ns %-13s PCE %.3f conv %d/%d res %d\n',tn,conds{icd},m.PCE,nOk,nV,nRes);
  end
  T = cell2table(rows,'VariableNames',{'tauN_ns','Condition','Jsc','Voc_mV','FF_pct','PCE_pct','nConverged','nResidualOK','nV'});
  writetable(T,fullfile(outDir,sprintf('ntype_fj_taun_part%d.csv',part)));
end
end

function [m,nOk,nRes,nV] = mirroredSweep(dev,G,Vgrid,Pin)
N = dev.N; Vg = Vgrid(:); nV = numel(Vg); idxStart = find(abs(Vg-0.02)<1e-12,1); S = cell(nV,1);
ws = @(s) struct('psiInit',s.psi,'phinInit',s.phin,'phipInit',s.phip);
r002 = gummelSolveQFPAnalytic(dev,-0.02,G,struct('psiInit',dev.psiEq,'phinInit',zeros(N,1),'phipInit',zeros(N,1)));
prev = r002; for v = [0.01,0.005,0], r0 = gummelSolveQFPAnalytic(dev,-v,G,ws(prev)); prev = r0; end
S{1} = r0; S{idxStart} = r002; prev = r002;
for k = idxStart+1:nV
  rk = gummelSolveQFPAnalytic(dev,-Vg(k),G,ws(prev)); S{k} = rk; if rk.isFinite, prev = rk; end
end
J = nan(nV,1); nOk = 0; nRes = 0;
for k = 1:nV
  r = S{k}; if isempty(r) || ~r.isFinite, continue; end
  J(k) = -r.J; nOk = nOk + (isfield(r,'numericalConverged') && r.numericalConverged);
  cr = checkContinuityResiduals(dev,G,r.n,r.p,r.JnSeg,r.JpSeg,r.psi); nRes = nRes + (cr.electron.pass && cr.hole.pass);
end
f = isfinite(J); m = pvMetrics([Vg(f), J(f)],Pin);
end
