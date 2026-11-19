function T = run_paper_r4_analyses(part, iSlice, nSlice)
%RUN_PAPER_R4_ANALYSES Analyses requested by the fourth simulated review (docs/HAKEM_SIMULASYONU_2026-11-07.md, roadmap group B).
% part 'L200': design map and rear-contact controls at CIE 6500 K / 200 lx (R3-1): for tau_p = 1-600 ns, RBSF-GaP and RHJ (r = 1-60) at
%              S_back = 1e5 cm/s, and HJ-0, HJ-U, F-HJ, RBSF-GaP at S_back = 1e5 and 1e3 cm/s (r = 1).
% part 'OPT' : per-architecture thickness/doping grid in the nanosecond regime (R2-7, R5-C3): tau = 1 and 10 ns (RHJ r = 1), S_back = 1e5,
%              absorber thickness scale 0.25-1.5 (nominal 1.5 um) x absorber doping 1e16-1e18 cm-3, HJ-0, HJ-U, F-HJ, RBSF-GaP, RHJ, AM1.5G and 1000 lx.
% part 'SRH' : injection-dependent SRH (R2-4): the SAME single SRH centre (tau_n0, tau_p0, E_t - E_i) in the p-type (RBSF-GaP) and the
%              n-type (RHJ) absorber instead of a constant lifetime ratio; tau pairs (ns) (10,10) (10,100) (100,10) (1,10) (10,1), E_t - E_i =
%              -0.3, 0, +0.3 eV; AM1.5G, 1000 and 200 lx. The design-map ratio is r = tau_p0/tau_n0 at tau_p = tau_n0.
% part 'SRH2': ratio scan r = tau_p0/tau_n0 = 2-30 at tau_n0 = 10 and 100 ns (E_t = E_i) and shallow levels E_t - E_i = +-0.8 eV at r = 10.
% Production physics (Stage-A ARC, SRH + radiative + Auger, mu_n = 500, chi_GaP = 3.8), 0.02 V grid, unified extractor (all finite points).
% Jobs are split into nSlice interleaved slices; output results/paper_r4_analyses/<part>_<iSlice>of<nSlice>.csv (checkpointed).
if nargin < 2, iSlice = 1; nSlice = 1; end
assertCanonicalPath(); warning('off','all');
root = fileparts(fileparts(mfilename('fullpath'))); outDir = fullfile(root,'results','paper_r4_analyses'); if ~isfolder(outDir), mkdir(outDir); end
arc = load(fullfile(root,'results','paper_stageA','stageA_comparison.mat'),'d1nm','d2nm');
arcFn = @(lam) tmmReflectance(lam,[1.0,1.38,2.35,gapRefractiveIndexReal(lam)],[arc.d1nm,arc.d2nm]);
Bmat = paperRadiativeCoefficients(300);
condNames = {'AM1.5G','6500K_1000lx','6500K_200lx'};
specs = {loadRealAM15G(), buildCIEIndoorSpectrum('6500K',1000), buildCIEIndoorSpectrum('6500K',200)};
Vgrids = {(0:0.02:1.70)', (0:0.02:1.56)', (0:0.02:1.50)'};
Pins = cellfun(@spectrumPower,specs);
% job list: {code, tauP_ns, ratio, Sback, tScale, Nabs, icond}
J = {};
switch part
  case 'L200'
    for tp = [1 3 10 30 100 300 600]
      J(end+1,:) = {'RBSF-GaP',tp,1,1e5,1,1e17,3}; %#ok<AGROW>
      for r = [1 1.5 2 3 5 10 20 30 60], J(end+1,:) = {'RHJ',tp,r,1e5,1,1e17,3}; end %#ok<AGROW>
      for c = {'HJ-0','HJ-U','F-HJ','RBSF-GaP'}, J(end+1,:) = {c{1},tp,1,1e3,1,1e17,3}; end %#ok<AGROW>
      for c = {'HJ-0','HJ-U','F-HJ'}, J(end+1,:) = {c{1},tp,1,1e5,1,1e17,3}; end %#ok<AGROW>
    end
  case 'OPT'
    for tp = [1 10]
      for c = {'HJ-0','HJ-U','F-HJ','RBSF-GaP','RHJ'}
        for ts = [0.25 0.5 0.75 1 1.5]
          for Na = [1e16 3e16 1e17 3e17 1e18]
            for ic = 1:2, J(end+1,:) = {c{1},tp,1,1e5,ts,Na,ic}; end %#ok<AGROW>
          end
        end
      end
    end
  case 'OPT2'   % extension of OPT beyond its doping/thickness edge (optimum at 1e18 cm-3 / 2.25 um)
    for tp = [1 10]
      for c = {'HJ-0','HJ-U','F-HJ','RBSF-GaP','RHJ'}
        for ts = [0.5 0.75 1 1.5 2]
          for Na = [3e18 1e19]
            for ic = 1:2, J(end+1,:) = {c{1},tp,1,1e5,ts,Na,ic}; end %#ok<AGROW>
          end
          if ts == 2, for Na = [1e17 3e17 1e18], for ic = 1:2, J(end+1,:) = {c{1},tp,1,1e5,ts,Na,ic}; end, end, end %#ok<AGROW>
        end
      end
    end
  case 'SRH'
    for pr = [10 10; 10 100; 100 10; 1 10; 10 1]'
      for Et = [-0.3 0 0.3]
        for c = {'RBSF-GaP','RHJ'}
          for ic = 1:3, J(end+1,:) = {c{1},pr(1),pr(2)/pr(1),1e5,1,1e17,ic,[pr(1) pr(2) Et]}; end %#ok<AGROW>
        end
      end
    end
  case 'SRH2'   % ratio scan of the single-centre case (crossover with the same centre in both absorbers) and shallow levels (+-0.8 eV)
    for tn = [10 100]
      for rr = [2 3 5 10 20 30]
        for c = {'RBSF-GaP','RHJ'}, for ic = 1:3, J(end+1,:) = {c{1},tn,rr,1e5,1,1e17,ic,[tn tn*rr 0]}; end, end %#ok<AGROW>
      end
      for Et = [-0.8 0.8]
        for c = {'RBSF-GaP','RHJ'}, for ic = 1:3, J(end+1,:) = {c{1},tn,10,1e5,1,1e17,ic,[tn tn*10 Et]}; end, end %#ok<AGROW>
      end
    end
  case 'DIM'    % round-5 review (R5-M2): design map at absorber thickness scale 0.5 and 1.5 (W = 0.75 and 2.25 um), mu_n = 500,
                % for a thickness-normalised (L/W) crossover criterion; RBSF-GaP (tau = tau_p) and RHJ (tau = r*tau_p), S_back = 1e5
    for ts = [0.5 1.5]
      for tp = [1 3 10 30 100 600]
        for ic = 1:2
          J(end+1,:) = {'RBSF-GaP',tp,1,1e5,ts,1e17,ic}; %#ok<AGROW>
          for r = [0.5 0.7 1 1.5 2 3 5 10 20 30 60 100], J(end+1,:) = {'RHJ',tp,r,1e5,ts,1e17,ic}; end %#ok<AGROW>
        end
      end
    end
  case 'BERN'   % re-run of a subset of OPT (nominal thickness/doping) and L200 jobs after the expm1 fix of sgBernoulli (round-5 review)
    for tp = [1 10]
      for c = {'HJ-0','HJ-U','F-HJ','RBSF-GaP','RHJ'}, for ic = 1:2, J(end+1,:) = {c{1},tp,1,1e5,1,1e17,ic}; end, end %#ok<AGROW>
    end
    for tp = [10 600]
      J(end+1,:) = {'RBSF-GaP',tp,1,1e5,1,1e17,3}; %#ok<AGROW>
      for r = [1 10], J(end+1,:) = {'RHJ',tp,r,1e5,1,1e17,3}; end %#ok<AGROW>
      for c = {'HJ-0','HJ-U','F-HJ'}, J(end+1,:) = {c{1},tp,1,1e5,1,1e17,3}; end %#ok<AGROW>
    end
  otherwise
    error('unknown part %s', part);
end
if size(J,2) < 8, J(:,8) = {[]}; end
J = J(iSlice:nSlice:end,:);
ck = fullfile(outDir,sprintf('%s_%dof%d.csv',part,iSlice,nSlice)); rows = {};
for ij = 1:size(J,1)
  [code,tp,r,Sb,ts,Na,ic,srh] = J{ij,:}; tauAbs = tp*1e-9*r; t0 = tic; Et = NaN;
  if strcmp(code,'RBSF-GaP') && isempty(srh), tauAbs = tp*1e-9; end
  ctr = [];
  if ~isempty(srh), ctr = struct('tauN',srh(1)*1e-9,'tauP',srh(2)*1e-9,'EtRelEi_eV',srh(3)); Et = srh(3); tauAbs = NaN; end
  try
    dev = buildPaperArchitecture(code,'AbsorberTau_s',max(tauAbs,1e-12),'AbsorberSRHCenters',ctr,'Sback',Sb,'AbsorberThicknessScale',ts,'AbsorberDoping',Na, ...
      'BradAbs',Bmat.GaAsPN_cm3ps,'BradGaP',Bmat.GaP_cm3ps,'BradGaNP',Bmat.GaPN_Ege18_cm3ps, ...
      'CnAugerAbs',3e-30,'CpAugerAbs',3e-30,'CnAugerGaP',3e-30,'CpAugerGaP',3e-30,'CnAugerGaNP',3e-30,'CpAugerGaNP',3e-30);
    G = deviceGenAndTransmitML(dev,specs{ic},arcFn);
    if ic == 1
      [~,~,m] = runAcceptedJVSweep(dev,G,Vgrids{ic},Pins(ic),true); cls = 'AM15G';
      Jsc = m.Jsc*1e3; Voc = m.Voc*1e3; FF = m.FF*100; PCE = m.PCE;
    else
      gate = indoorSweepLocal(dev,G,Vgrids{ic},Pins(ic)); cls = char(gate.curveClass);
      Jsc = gate.Jsc*1e6; Voc = gate.Voc*1e3; FF = gate.FF*100; PCE = gate.PCE;
    end
  catch ME
    Jsc = NaN; Voc = NaN; FF = NaN; PCE = NaN; cls = ['ERR ' strrep(ME.message,',',';')];
  end
  rows(end+1,:) = {code, tp, r, Sb, ts, Na, condNames{ic}, tauAbs, Et, Jsc, Voc, FF, PCE, cls, toc(t0)}; %#ok<AGROW>
  T = cell2table(rows,'VariableNames',{'Architecture','tauP_ns','ratio','Sback','tScale','Nabs','Condition','tauAbs_s','EtRelEi_eV','Jsc','Voc_mV','FF_pct','PCE_pct','Class','Seconds'});
  writetable(T,ck);
  fprintf('%d/%d %s tau=%g r=%g Sb=%g t=%g N=%g %s PCE=%.3f (%.0f s)\n',ij,size(J,1),code,tp,r,Sb,ts,Na,condNames{ic},PCE,toc(t0));
end
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
fU = gate.isFinite(:); mU = pvMetrics([gate.Vs(fU), gate.Jterm(fU)],Pin);   % unified extractor: all finite points, J_sc at 0 V, linear V_oc, parabolic MPP
gate.Jsc = mU.Jsc; gate.Voc = mU.Voc; gate.FF = mU.FF; gate.PCE = mU.PCE;
end
