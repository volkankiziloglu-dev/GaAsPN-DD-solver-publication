function T = run_ntype_front_junction_check(tauList_ns, tag)
%RUN_NTYPE_FRONT_JUNCTION_CHECK Does the RHJ gain come from the n-type absorber or from the rear junction?
% Adds the missing combination "n-type absorber, junction at the front" (RBSF-GaP-n: p+-GaP window / n-GaAsPN / n+-GaP rear BSF)
% to the lifetime table of the paper (tab:lifetime): p-type lifetime tau_p, n-type lifetime r*tau_p for RHJ and RBSF-GaP-n.
% RBSF-GaP (tau_p) and RHJ (r*tau_p) are recomputed as a reproduction check. Settings, sweeps and extractors as in
% run_paper_design_map.m (unified extractor, S_back = 1e5 cm/s, AM1.5G and CIE 6500 K / 1000 lx).
% Output: results/paper_ntype_front_junction/ntype_fj_<tag>.csv
if nargin < 1 || isempty(tauList_ns), tauList_ns = [600 100 10 1]; end
if nargin < 2, tag = 'all'; end
assertCanonicalPath(); warning('off','all');
root = fileparts(fileparts(mfilename('fullpath'))); outDir = fullfile(root,'results','paper_ntype_front_junction'); if ~isfolder(outDir), mkdir(outDir); end
arc = load(fullfile(root,'results','paper_stageA','stageA_comparison.mat'),'d1nm','d2nm');
arcFn = @(lam) tmmReflectance(lam,[1.0,1.38,2.35,gapRefractiveIndexReal(lam)],[arc.d1nm,arc.d2nm]);
Bmat = paperRadiativeCoefficients(300);
conds = {'AM1.5G'; '6500K_1000lx'}; specs = {loadRealAM15G(); buildCIEIndoorSpectrum('6500K',1000)};
Vgrids = {(0:0.02:1.70)'; (0:0.02:1.56)'}; Pins = cellfun(@spectrumPower,specs);
rList = [1 3 10 100]; rows = {}; ck = fullfile(outDir,sprintf('ntype_fj_%s.csv',tag));
for tp = tauList_ns(:)'
  jobs = {'RBSF-GaP',1};
  for r = rList, jobs(end+1,:) = {'RHJ',r}; jobs(end+1,:) = {'RBSF-GaP-n',r}; end %#ok<AGROW>
  for icd = 1:2
    spec = specs{icd}; Vgrid = Vgrids{icd}; Pin = Pins(icd);
    for ij = 1:size(jobs,1)
      code = jobs{ij,1}; r = jobs{ij,2}; tauAbs = tp*1e-9*r; t0 = tic;
      try
        dev = buildPaperArchitecture(code,'AbsorberTau_s',tauAbs,'Sback',1e5, ...
          'BradAbs',Bmat.GaAsPN_cm3ps,'BradGaP',Bmat.GaP_cm3ps,'BradGaNP',Bmat.GaPN_Ege18_cm3ps, ...
          'CnAugerAbs',3e-30,'CpAugerAbs',3e-30,'CnAugerGaP',3e-30,'CpAugerGaP',3e-30,'CnAugerGaNP',3e-30,'CpAugerGaNP',3e-30);
        G = deviceGenAndTransmitML(dev,spec,arcFn);
        if code == "RBSF-GaP-n"
          % p-type front: forward bias is V < 0 and the photocurrent has the opposite sign (tested). Mirrored sweep
          % (V -> -V, J -> -J) with the same solver, continuation and unified extractor; the acceptance gate assumes the n-front sign
          % convention, so each point is checked with the native continuity residuals instead.
          [m,nOk,nRes] = mirroredSweep(dev,G,Vgrid,Pin); cls = sprintf('mirror conv %d/%d resOK %d',nOk,numel(Vgrid),nRes);
          Jsc = m.Jsc*1e3*(icd==1) + m.Jsc*1e6*(icd==2); Voc = m.Voc*1e3; FF = m.FF*100; PCE = m.PCE;
        elseif icd == 1
          [~,~,m] = runAcceptedJVSweep(dev,G,Vgrid,Pin,true); cls = 'AM15G';
          Jsc = m.Jsc*1e3; Voc = m.Voc*1e3; FF = m.FF*100; PCE = m.PCE;
        else
          gate = indoorSweepLocal(dev,G,Vgrid,Pin); cls = char(gate.curveClass);
          Jsc = gate.Jsc*1e6; Voc = gate.Voc*1e3; FF = gate.FF*100; PCE = gate.PCE;
        end
      catch ME
        Jsc = NaN; Voc = NaN; FF = NaN; PCE = NaN; cls = ['ERR ' strrep(ME.message,',',';')];
      end
      rows(end+1,:) = {tp, r, conds{icd}, code, tauAbs, Jsc, Voc, FF, PCE, cls}; %#ok<AGROW>
      fprintf('%5g ns r=%-4g %-13s %-11s Jsc %.3f Voc %.1f FF %.2f PCE %.3f %s (%.0f s)\n',tp,r,conds{icd},code,Jsc,Voc,FF,PCE,cls,toc(t0));
    end
  end
  T = cell2table(rows,'VariableNames',{'tauP_ns','ratio','Condition','Architecture','tauAbs_s','Jsc','Voc_mV','FF_pct','PCE_pct','Class'});
  writetable(T,ck);
end
end

function [m,nOk,nRes] = mirroredSweep(dev,G,Vgrid,Pin)
N = dev.N; Vg = Vgrid(:); idxStart = find(abs(Vg-0.02)<1e-12,1); S = cell(numel(Vg),1);
ws = @(s) struct('psiInit',s.psi,'phinInit',s.phin,'phipInit',s.phip);
r002 = gummelSolveQFPAnalytic(dev,-0.02,G,struct('psiInit',dev.psiEq,'phinInit',zeros(N,1),'phipInit',zeros(N,1)));
prev = r002; for v = [0.01,0.005,0], r0 = gummelSolveQFPAnalytic(dev,-v,G,ws(prev)); prev = r0; end
S{1} = r0; S{idxStart} = r002; prev = r002;
for k = idxStart+1:numel(Vg)
  rk = gummelSolveQFPAnalytic(dev,-Vg(k),G,ws(prev)); S{k} = rk; if rk.isFinite, prev = rk; end
end
J = nan(numel(Vg),1); nOk = 0; nRes = 0;
for k = 1:numel(Vg)
  r = S{k}; if isempty(r) || ~r.isFinite, continue; end
  J(k) = -r.J; nOk = nOk + (isfield(r,'numericalConverged') && r.numericalConverged);
  cr = checkContinuityResiduals(dev,G,r.n,r.p,r.JnSeg,r.JpSeg,r.psi); nRes = nRes + (cr.electron.pass && cr.hole.pass);
end
f = isfinite(J); m = pvMetrics([Vg(f), J(f)],Pin);
end

function gate = indoorSweepLocal(dev,G,Vgrid,Pin)
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
