function T = run_lowlight_experiment(method, codes, luxList, opts)
%RUN_LOWLIGHT_EXPERIMENT EXPERIMENTAL: two remedies for the low-light (CIE 6500 K, 100/200 lx) SCREEN problem.
% method:
%   'prod'          production solver (gummelSolveQFPAnalytic) -- baseline
%   'expm1'         same solver, only the sgBernoulliExpm1 kernel in the flux
%   'excess'        excess-carrier (equilibrium-subtracted) SG flux (assembleResidualQFPX 'excess')
%   'prod+vpaeval'  'prod' states, fluxes/residuals RE-EVALUATED in vpa (no polishing)
%   'prod+vpa'      'prod' states + mixed-precision Newton polishing with vpa residual (production discretization)
%   'excess+dd'     'excess' states + excess-carrier polishing with hi/lo (double-double) state
% Device setup identical to production (run_paper_voc_lux.m / run_paper_spatial_profiles.m): default
% AbsorberTau 600 ns, BradAbs/GaP/GaNP = paperRadiativeCoefficients(300), Auger 3e-30, ARC stageA,
% buildCIEIndoorSpectrum('6500K',lux), Vgrid 0:0.02:1.56, V=0 via reverse continuation, PRODUCTION acceptance
% gate indoorAcceptanceGate.m (unchanged). Outputs: results/experimental_lowlight/.
if nargin < 2 || isempty(codes), codes = ["RBSF-GaP","RHJ"]; end
if nargin < 3 || isempty(luxList), luxList = [100 200]; end
if nargin < 4, opts = struct(); end
if ~isfield(opts,'nDigits'), opts.nDigits = 40; end
if ~isfield(opts,'maxSweeps'), opts.maxSweeps = 6; end
if ~isfield(opts,'convFlag'), opts.convFlag = 'polish'; end  % convergence flag for the polished state: 'polish' (own Newton step) | 'inherit' (from the base solution)
if ~isfield(opts,'saveStates'), opts.saveStates = true; end
if ~isfield(opts,'spectrum'), opts.spectrum = '6500K'; end    % '6500K' (CIE, luxList) | 'AM15G' (luxList ignored)
isAM = strcmpi(opts.spectrum,'AM15G');
if isAM, luxList = NaN; end
if ~isfield(opts,'meshScale'), opts.meshScale = 1; end   % buildPaperArchitecture meshScale (mesh regression)
ms = opts.meshScale; mtag = ''; if ms ~= 1, mtag = sprintf('_m%g',ms); end
warning('off','all');
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root,'src','experimental_lowlight'));
outDir = fullfile(root,'results','experimental_lowlight'); stDir = fullfile(outDir,'states');
if ~isfolder(stDir), mkdir(stDir); end
arc = load(fullfile(root,'results','paper_stageA','stageA_comparison.mat'),'d1nm','d2nm');
arcFn = @(lam) tmmReflectance(lam,[1.0,1.38,2.35,gapRefractiveIndexReal(lam)],[arc.d1nm,arc.d2nm]);
Bmat = paperRadiativeCoefficients(300);
prodTab = readtable(fullfile(root,'results','paper_voc_lux','lux_series_pce.csv'),'TextType','string');
Vgrid = (0:0.02:1.56)';
if isAM
  % AM1.5G reference: run_paper_absorber_doping_ranking (AM_UNIFIED=1), 1e17 rows = paper Table (outdoor)
  Vgrid = (0:0.02:1.70)';
  amRef = readtable(fullfile(root,'results','paper_absorber_doping_unified_am','doping_ranking.csv'),'TextType','string');
  amRef = amRef(abs(amRef.N_cm3-1e17) < 1, :);
  prodTab = table(amRef.Architecture, nan(height(amRef),1), amRef.Jsc*1e3, amRef.Voc_mV, amRef.FF_pct, amRef.PCE_pct, amRef.Class, ...
    'VariableNames',{'Architecture','lux','Jsc_uAcm2','Voc_mV','FF_pct','PCE_pct','Class'});
end
condSuffix = '';
if ms ~= 1   % mesh reference: results/paper_mesh_richardson (production, 1000 lx and AM1.5G; x2/x4 columns)
  rs = readtable(fullfile(root,'results','paper_mesh_richardson','richardson_summary.csv'),'TextType','string');
  if isAM, rs = rs(rs.Condition=="AM1.5G",:); else, rs = rs(rs.Condition=="6500K_1000lx",:); end
  col = sprintf('x%g',ms); archs = unique(rs.Architecture); refRows = {};
  for a = archs'
    g = @(mm) rs.(col)(rs.Architecture==a & rs.metric==mm);
    mc = readtable(fullfile(root,'results','paper_mesh_richardson',sprintf('mesh_richardson_%s.csv',a)),'TextType','string');
    cl = mc.Class(mc.MeshScale==ms & contains(mc.Condition,ifelse_(isAM,'AM1.5G','1000lx')));
    refRows(end+1,:) = {a, 1000, g("Jsc")*1e3, g("Voc"), g("FF"), g("PCE"), cl(1)}; %#ok<AGROW>
  end
  prodTab = cell2table(refRows,'VariableNames',{'Architecture','lux','Jsc_uAcm2','Voc_mV','FF_pct','PCE_pct','Class'});
end
if isAM, condSuffix = '__AM15G';
elseif ~isequal(sort(luxList(:))',[100 200]), condSuffix = ['__' strjoin(compose('%glx',luxList(:)'),'_')]; end
condSuffix = [condSuffix mtag];
tag = strrep(strrep(method,'+','_'),'>','_to_');
if contains(method,'+') && strcmp(opts.convFlag,'inherit'), tag = [tag '_miras']; end
if isfield(opts,'stepTol'), tag = [tag sprintf('_tol%g',opts.stepTol)]; end   % polishing flag threshold (polishStateHPv3 opts.stepTol)
rows = {};
for code = string(codes(:))'
  dev = buildPaperArchitecture(code, meshScale=ms, BradAbs=Bmat.GaAsPN_cm3ps, BradGaP=Bmat.GaP_cm3ps, BradGaNP=Bmat.GaPN_Ege18_cm3ps, ...
    CnAugerAbs=3e-30,CpAugerAbs=3e-30,CnAugerGaP=3e-30,CpAugerGaP=3e-30,CnAugerGaNP=3e-30,CpAugerGaNP=3e-30);
  for lux = luxList(:)'
    t0 = tic;
    if isAM, spec = loadRealAM15G(); cond = 'AM15G';
    else, spec = buildCIEIndoorSpectrum('6500K',lux); cond = sprintf('%glx',lux); end
    cond = [cond mtag];
    Pin = spectrumPower(spec);
    G = deviceGenAndTransmitML(dev,spec,arcFn);
    stFile = @(m) fullfile(stDir, sprintf('%s_%s_%s.mat', code, cond, strrep(m,'>','_to_')));
    polishInfo = [];
    switch method
      case 'prod'
        if isAM, states = sweepAM(@(V,init) gummelSolveQFPAnalytic(dev,V,G,init), dev, Vgrid);
        else, states = sweepLocal(@(V,init) gummelSolveQFPAnalytic(dev,V,G,init), dev, Vgrid); end
      case {'expm1','excess','excessB0'}
        devX = dev; devX.fluxModeX = method;
        states = sweepLocal(@(V,init) gummelSolveQFPAnalyticX(devX,V,G,init), dev, Vgrid);
      case {'excess_rs','excessB0_rs'}
        devX = dev; devX.fluxModeX = extractBefore(method,'_rs');
        states = sweepLocal(@(V,init) gummelSolveQFPAnalyticX(devX,V,G,setfield(init,'rowScale',true)), dev, Vgrid); %#ok<SFLD>
      case {'prod>excess','prod>excessB0','prod>expm1'}
        % warm start: at each V, start from the converged production state and solve with the excess-carrier residual
        devX = dev; devX.fluxModeX = extractAfter(method,'>');
        L = load(stFile('prod'),'states'); P0 = L.states; states = P0;
        for k = 1:numel(P0)
          rk = gummelSolveQFPAnalyticX(devX,Vgrid(k),G,struct('psiInit',P0(k).psi,'phinInit',P0(k).phin,'phipInit',P0(k).phip));
          states(k) = rk;
        end
      case {'prod+vpaeval','prod+vpa','excess+dd','prod>excess+dd','prod+dd','prod+dd3'}
        base = extractBefore(method,'+'); L = load(stFile(base),'states'); states = L.states;
        switch method
          case 'prod+vpaeval', kind = 'vpadiff'; ns = 0;
          case 'prod+vpa',     kind = 'vpadiff'; ns = opts.maxSweeps;
          case {'excess+dd','prod>excess+dd','prod+dd','prod+dd3'}, kind = 'ddexcess'; ns = opts.maxSweeps;   % prod+dd: directly from the production state, without the excess warm start
        end
        polishInfo = cell(numel(states),1);
        for k = 1:numel(states)
          if ~states(k).isFinite, continue; end
          popts = struct('maxSweeps',ns,'nDigits',opts.nDigits,'convFlag',opts.convFlag);
          if isfield(opts,'stepTol'), popts.stepTol = opts.stepTol; end
          if strcmp(method,'prod+dd3'), [sk, polishInfo{k}] = polishStateHPv3(dev, states(k), Vgrid(k), G, kind, popts);
          else, [sk, polishInfo{k}] = polishStateHPv2(dev, states(k), Vgrid(k), G, kind, popts); end
          states(k) = orderfields(sk, states(k));
        end
      otherwise
        error('unknown method %s', method);
    end
    if opts.saveStates && any(strcmp(method,{'prod','expm1','excess','excessB0','prod>excess','prod>expm1'})), save(stFile(method),'states','Vgrid','-v7.3'); end
    gate = indoorAcceptanceGate(dev,G,Vgrid,states,Pin);
    % --- point-wise diagnostic (region-scaled rJ: same as rJatK in the gate) ---
    rJk = nan(numel(Vgrid),1);
    for k = 1:numel(Vgrid)
      if ~gate.isFinite(k), continue; end
      switch gate.region(k)
        case {"shortcircuit","voc"}, s = abs(gate.JscRaw);
        case "mpp", s = max(abs(gate.Jterm(k)),0.1*abs(gate.JscRaw));
        otherwise, s = max(abs(gate.Jterm(k)),1e-9);
      end
      rJk(k) = (max(states(k).Jseg)-min(states(k).Jseg))/max(s,1e-9);
    end
    P = table(Vgrid,gate.region,gate.Jterm,rJk,gate.rFB,gate.rnGlobal,gate.rpGlobal, ...
      gate.numericalConverged,gate.poissonConverged,gate.pointPass, ...
      'VariableNames',{'V_V','region','Jterm_Acm2','rJ','rFB','rn','rp','numConv','poissonConv','pass'});
    writetable(P, fullfile(outDir, sprintf('points_%s_%s_%s.csv', code, cond, tag)));
    f = gate.isFinite(:);
    writetable(table(Vgrid(f),gate.Jterm(f),double(gate.pointPass(f)),'VariableNames',{'V_V','J_Acm2','Accepted'}), ...
      fullfile(outDir, sprintf('jv_%s_%s_%s.csv', code, cond, tag)));
    crit = f & gate.region ~= "deepforward";
    if isAM, m = pvMetrics([Vgrid(f), gate.Jterm(f)], Pin); else, m = pvMetricsUnified(Vgrid(f), gate.Jterm(f), Pin); end
    nQfp = sum(arrayfun(@(s) qfpResultAccepted(s), states(f)));
    if isAM, pr = prodTab(prodTab.Architecture==code,:); else, pr = prodTab(prodTab.Architecture==code & prodTab.lux==lux,:); end
    if isempty(pr), pr = table(NaN,NaN,NaN,NaN,"",'VariableNames',{'Jsc_uAcm2','Voc_mV','FF_pct','PCE_pct','Class'}); end
    fails = @(v,th) sum(crit & ~(v < th));
    rows(end+1,:) = {char(code), lux, cond, method, char(gate.curveClass), gate.scPass, gate.vocPass, gate.mppPass, ...
      gate.leftPass, gate.rightPass, sum(gate.pointPass(crit)), sum(crit), ...
      max(rJk(crit)), max(gate.rFB(crit)), max(gate.rnGlobal(crit)), max(gate.rpGlobal(crit)), ...
      median(rJk(crit)), median(gate.rFB(crit)), median(gate.rnGlobal(crit)), median(gate.rpGlobal(crit)), ...
      fails(rJk,0.010), fails(gate.rFB,0.005), fails(gate.rnGlobal,0.15), fails(gate.rpGlobal,0.015), ...
      sum(crit & ~(gate.numericalConverged(:) & gate.poissonConverged(:))), ...
      gate.Jnoise, gate.JphotoScale, m.Jsc*1e6, m.Voc*1e3, m.FF*100, m.PCE, ...
      m.Jsc*1e6-pr.Jsc_uAcm2, m.Voc*1e3-pr.Voc_mV, m.FF*100-pr.FF_pct, m.PCE-pr.PCE_pct, char(pr.Class), nQfp, sum(f), toc(t0)}; %#ok<AGROW>
    fprintf('%s %s %-13s : %s  (critical passed %d/%d)  max rJ %.3g rFB %.3g rn %.3g rp %.3g | %.1f s\n', code, cond, method, ...
      gate.curveClass, sum(gate.pointPass(crit)), sum(crit), max(rJk(crit)), max(gate.rFB(crit)), ...
      max(gate.rnGlobal(crit)), max(gate.rpGlobal(crit)), toc(t0));
    if strcmp(method,'prod+dd3')
      ok = find(~cellfun(@isempty,polishInfo)); PI = [polishInfo{ok}];
      D = table(Vgrid(ok), [PI.v2Converged]', [PI.polishConverged]', [PI.startOK]', [PI.totalCorr]', [PI.certStep]', [PI.bestIdx]', [PI.sweeps]', gate.numericalConverged(ok), ...
        'VariableNames',{'V_V','v2Flag','certFlag','startOK','totalCorr_V','certStep_V','bestIdx','sweeps','numConvFinal'});
      writetable(D, fullfile(outDir, sprintf('polishdiag_%s_%s_%s.csv', code, cond, tag)));
      fprintf('   v3 flag: converged %d/%d (v2 %d), startOK %d, max total correction %.3g V, max certStep %.3g V\n', ...
        sum(D.numConvFinal), height(D), sum(D.v2Flag), sum(D.startOK), max(D.totalCorr_V), max(D.certStep_V));
    end
    if ~isempty(polishInfo)
      hh = cellfun(@(q) q.normHistory(end,:)./max(q.normHistory(1,:),realmin), polishInfo(~cellfun(@isempty,polishInfo)), 'UniformOutput', false);
      hh = vertcat(hh{:}); fprintf('   polishing: |F_n| last/first median %.3g, |F_p| last/first median %.3g\n', median(hh(:,1)), median(hh(:,2)));
    end
  end
end
T = cell2table(rows,'VariableNames',{'Architecture','lux','Condition','Method','Class','scPass','vocPass','mppPass','leftPass','rightPass', ...
  'nPassCrit','nCrit','max_rJ','max_rFB','max_rn','max_rp','med_rJ','med_rFB','med_rn','med_rp', ...
  'nFail_rJ','nFail_rFB','nFail_rn','nFail_rp','nFail_conv','Jnoise_Acm2','JphotoScale_Acm2', ...
  'Jsc_uAcm2','Voc_mV','FF_pct','PCE_pct','dJsc_uAcm2','dVoc_mV','dFF_pct','dPCE_pct','ProdClass','nQfpAccepted','nFinite','Seconds'});
writetable(T, fullfile(outDir, sprintf('summary_%s__%s%s.csv', tag, strjoin(string(codes),'_'), condSuffix)));
end

function states = sweepLocal(solve, dev, Vgrid)
% Same continuation scheme as indoorSweepLocal in run_paper_voc_lux.m.
N = dev.N; idxStart = find(abs(Vgrid-0.02)<1e-12,1);
r002 = solve(0.02, struct('psiInit',dev.psiEq,'phinInit',zeros(N,1),'phipInit',zeros(N,1)));
prev = r002;
for v = [0.01,0.005,0], r0 = solve(v, struct('psiInit',prev.psi,'phinInit',prev.phin,'phipInit',prev.phip)); prev = r0; end
states = repmat(r002,numel(Vgrid),1); states(1) = r0; states(idxStart) = r002; prev = r002;
for k = idxStart+1:numel(Vgrid)
  rk = solve(Vgrid(k), struct('psiInit',prev.psi,'phinInit',prev.phin,'phipInit',prev.phip));
  states(k) = rk; if rk.isFinite, prev = rk; end
end
end

function m = pvMetricsUnified(V, J, Pin)
% MATLAB counterpart of pv_metrics in cases/unified_metrics.py (extractor of the production lux table).
[V,o] = sort(V(:)); J = J(o); J = J(:);
Jsc = NaN;
for k = 1:numel(V)-1
  if V(k) <= 0 && 0 <= V(k+1) && V(k+1) > V(k), t = (0-V(k))/(V(k+1)-V(k)); Jsc = J(k)+t*(J(k+1)-J(k)); break; end
  if V(k) == 0, Jsc = J(k); break; end
end
Voc = NaN;
for k = 1:numel(V)-1
  if V(k) >= 0 && J(k) >= 0 && J(k+1) < 0, t = (0-J(k))/(J(k+1)-J(k)); Voc = V(k)+t*(V(k+1)-V(k)); break; end
end
if ~(isfinite(Jsc) && isfinite(Voc) && Voc > 0 && Jsc > 0), m = struct('Jsc',Jsc,'Voc',Voc,'FF',NaN,'PCE',NaN); return; end
Pmax = 0; kb = 0;
for k = 1:numel(V)
  if V(k) < 0 || V(k) > Voc || J(k) < 0, continue; end
  if V(k)*J(k) > Pmax, Pmax = V(k)*J(k); kb = k; end
end
if kb > 1 && kb < numel(V)
  V1 = V(kb-1); V2 = V(kb); V3 = V(kb+1); P1 = V1*J(kb-1); P2 = Pmax; P3 = V3*J(kb+1);
  num = (V2-V1)^2*(P2-P3) - (V2-V3)^2*(P2-P1); den = (V2-V1)*(P2-P3) - (V2-V3)*(P2-P1);
  if abs(den) > 0
    Vv = V2 - 0.5*num/den;
    if V1 < Vv && Vv < V3
      Pv = P1*(Vv-V2)*(Vv-V3)/((V1-V2)*(V1-V3)) + P2*(Vv-V1)*(Vv-V3)/((V2-V1)*(V2-V3)) + P3*(Vv-V1)*(Vv-V2)/((V3-V1)*(V3-V2));
      if Pv > Pmax, Pmax = Pv; end
    end
  end
end
m = struct('Jsc',Jsc,'Voc',Voc,'FF',Pmax/(Voc*Jsc),'PCE',Pmax*1e4/Pin*100);
end

function states = sweepAM(solve, dev, Vgrid)
% Same continuation scheme as src/runAcceptedJVSweep.m (AM1.5G production sweep): starts
% at 0.02 V, reverse continuation to Vgrid<0.02 points, forward sweep; stops after a non-finite point.
N = dev.N; idxStart = find(abs(Vgrid-0.02)<1e-12,1);
r = solve(Vgrid(idxStart), struct('psiInit',dev.psiEq,'phinInit',zeros(N,1),'phipInit',zeros(N,1)));
states = repmat(r,numel(Vgrid),1); for k = 1:numel(Vgrid), states(k).isFinite = false; end
states(idxStart) = r;
if r.isFinite
  prev0 = r;
  for v = sort(Vgrid(Vgrid<0.02),'descend')'
    r0 = solve(v, struct('psiInit',prev0.psi,'phinInit',prev0.phin,'phipInit',prev0.phip));
    states(find(abs(Vgrid-v)<1e-12,1)) = r0; prev0 = r0;
  end
end
prev = r;
for k = idxStart+1:numel(Vgrid)
  if ~prev.isFinite, break; end
  rk = solve(Vgrid(k), struct('psiInit',prev.psi,'phinInit',prev.phin,'phipInit',prev.phip));
  states(k) = rk; if rk.isFinite, prev = rk; end
end
end

function v = ifelse_(c, a, b)
if c, v = a; else, v = b; end
end
