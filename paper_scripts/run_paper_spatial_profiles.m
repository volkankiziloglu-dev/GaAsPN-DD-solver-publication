function T = run_paper_spatial_profiles()
%RUN_PAPER_SPATIAL_PROFILES Electron, hole and total current density and carrier densities versus position at V = 0, at the maximum
% power point (grid maximum of V*J) and at open circuit, for the main-text architectures under AM1.5G and CIE 6500 K / 1000 lx.
% Follows the Sb2S3 cur_profile.m / np_profile.m pattern. Same devices, recombination, ARC and QFP continuation as the production runs
% (run_paper_stageA_comparison.m, run_paper_indoor_production.m). Current conservation: J_tot(x) = J_n + J_p over all segments is
% reported as max-min spread (the solver's own checkCurrentConservation metric). Outputs: results/paper_spatial_profiles/.
assertCanonicalPath();
warning('off','MATLAB:nearlySingularMatrix'); warning('off','MATLAB:singularMatrix');
root = fileparts(fileparts(mfilename('fullpath')));
outDir = fullfile(root,'results','paper_spatial_profiles'); if ~isfolder(outDir), mkdir(outDir); end
figDir = fullfile(root,'figures','main_results');
arc = load(fullfile(root,'results','paper_stageA','stageA_comparison.mat'),'d1nm','d2nm');
arcFn = @(lam) tmmReflectance(lam,[1.0,1.38,2.35,gapRefractiveIndexReal(lam)],[arc.d1nm,arc.d2nm]);
Bmat = paperRadiativeCoefficients(300);
codes = ["RBSF-GaP","HJ-0"];
conds = struct('name',{'AM1.5G','6500K_1000lx'}, 'label',{'AM1.5G','6500 K, 1000 lx'}, 'unit',{'mA cm^{-2}','\muA cm^{-2}'}, 'scale',{1e3,1e6});
Vgrid = (0:0.02:1.56)';
rows = {}; fid = fopen(fullfile(outDir,'spatial_profiles_report.txt'),'w'); c0 = onCleanup(@() fclose(fid));
for ic = 1:numel(codes)
  for icond = 1:numel(conds)
    cd_ = conds(icond);
    if icond == 1, spec = loadRealAM15G(); else, spec = buildCIEIndoorSpectrum('6500K',1000); end
    Pin = spectrumPower(spec);
    [dev,meta] = buildPaperArchitecture(codes(ic), BradAbs=Bmat.GaAsPN_cm3ps, BradGaP=Bmat.GaP_cm3ps, BradGaNP=Bmat.GaPN_Ege18_cm3ps, ...
      CnAugerAbs=3e-30,CpAugerAbs=3e-30,CnAugerGaP=3e-30,CpAugerGaP=3e-30,CnAugerGaNP=3e-30,CpAugerGaNP=3e-30);
    G = deviceGenAndTransmitML(dev,spec,arcFn);
    fprintf('\n%s | %s: sweep\n', codes(ic), cd_.label);
    [st, J] = continuationSweep(dev, G, Vgrid);
    ok = isfinite(J) & Vgrid > 0.02 & J > 0; [~,km] = max(Vgrid(ok).*J(ok)); idx = find(ok); km = idx(km);
    % open circuit: secant refinement between the last positive and first negative grid point
    kpos = find(J > 0 & isfinite(J), 1, 'last'); kneg = kpos + 1; V1 = Vgrid(kpos); J1 = J(kpos); V2 = Vgrid(kneg); J2 = J(kneg);
    prev = st{kpos}; rOc = prev; Vo = V1;
    for it = 1:6
      Vo = V1 - J1*(V2-V1)/(J2-J1);
      rOc = gummelSolveQFPAnalytic(dev,Vo,G,struct('psiInit',prev.psi,'phinInit',prev.phin,'phipInit',prev.phip));
      if abs(rOc.J) < 1e-4*abs(J(1)), break; end
      if rOc.J > 0, V1 = Vo; J1 = rOc.J; prev = rOc; else, V2 = Vo; J2 = rOc.J; end
    end
    sel = {st{1}, st{km}, rOc}; Vsel = [0, Vgrid(km), Vo];
    labs = {'V = 0 (J_{sc})', sprintf('V = %.2f V (MPP)',Vgrid(km)), sprintf('V = %.3f V (V_{oc})',Vo)};
    xnm = dev.x*1e7; xseg = 0.5*(xnm(1:end-1)+xnm(2:end)); lsn = dev.layerStartNode;
    f = figure('Visible','off','Units','centimeters','Position',[2 2 24 14],'Color','w');
    tl = tiledlayout(f,2,3,'TileSpacing','compact','Padding','compact');
    for k = 1:3
      s = sel{k}; Jn = s.JnSeg*cd_.scale; Jp = s.JpSeg*cd_.scale; Jt = Jn + Jp;
      sp = max(Jt) - min(Jt);
      ax = nexttile(tl,k); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
      for j = 2:numel(lsn), xline(ax,xnm(lsn(j))/1e3,':','Color',[.55 .55 .55],'HandleVisibility','off'); end
      yline(ax,0,'-','Color',[.8 .8 .8],'HandleVisibility','off');
      h1 = plot(ax,xseg/1e3,Jn,'-','Color',[0 .447 .741],'LineWidth',1.4);
      h2 = plot(ax,xseg/1e3,Jp,'-','Color',[.85 .325 .098],'LineWidth',1.4);
      h3 = plot(ax,xseg/1e3,Jt,'k--','LineWidth',1.1);
      title(ax,labs{k},'FontWeight','normal'); xlim(ax,[0 xnm(end)/1e3]);
      if k == 1, ylabel(ax,sprintf('J (%s)',cd_.unit)); legend(ax,[h1 h2 h3],{'J_n','J_p','J_{tot}'},'Location','east','Box','off'); end
      ax = nexttile(tl,k+3); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
      for j = 2:numel(lsn), xline(ax,xnm(lsn(j))/1e3,':','Color',[.55 .55 .55],'HandleVisibility','off'); end
      g1 = plot(ax,xnm/1e3,s.n,'-','Color',[0 .447 .741],'LineWidth',1.4); g2 = plot(ax,xnm/1e3,s.p,'-','Color',[.85 .325 .098],'LineWidth',1.4);
      set(ax,'YScale','log'); xlim(ax,[0 xnm(end)/1e3]); ylim(ax,[1e0 1e20]); yticks(ax,10.^(0:4:20)); xlabel(ax,'x (\mum)');
      if k == 1, ylabel(ax,'n, p (cm^{-3})'); legend(ax,[g1 g2],{'n','p'},'Location','southeast','Box','off'); end
      cc = s.currentConservation;
      fprintf(fid,'%s | %s | %s: Jterminal %.6g %s | Jtot spread %.3e %s (%.2e of |J|) | solver rJ %.2e | continuity e %.1e h %.1e\n', ...
        codes(ic), cd_.label, labs{k}, s.J*cd_.scale, cd_.unit, sp, cd_.unit, sp/max(abs(s.J*cd_.scale),1e-12), cc.rJ, ...
        s.continuityResiduals.electron.maxRelGlobal, s.continuityResiduals.hole.maxRelGlobal);
      rows(end+1,:) = {char(codes(ic)), cd_.name, Vsel(k), s.J*1e6, sp/cd_.scale*1e6, cc.rJ, ...
        s.continuityResiduals.electron.maxRelGlobal, s.continuityResiduals.hole.maxRelGlobal, min(s.n), min(s.p), max(s.n), max(s.p)}; %#ok<AGROW>
      vtags = {'V0','MPP','Voc'};
      tag = sprintf('%s_%s_%s',codes(ic),cd_.name,vtags{k});
      writetable(table(xnm(:),s.n(:),s.p(:),'VariableNames',{'x_nm','n_cm3','p_cm3'}),fullfile(outDir,[tag '_carriers.csv']));
      writetable(table(xseg(:),s.JnSeg(:)*1e6,s.JpSeg(:)*1e6,(s.JnSeg(:)+s.JpSeg(:))*1e6,'VariableNames',{'x_nm','Jn_uAcm2','Jp_uAcm2','Jtot_uAcm2'}),fullfile(outDir,[tag '_currents.csv']));
    end
    title(tl, sprintf('%s, %s',codes(ic),cd_.label),'FontWeight','normal');
    stem = sprintf('fig_spatial_%s_%s',strrep(char(codes(ic)),'-','_'),cd_.name);
    exportgraphics(f,fullfile(figDir,[stem '.png']),'Resolution',300); exportgraphics(f,fullfile(figDir,[stem '.pdf']),'ContentType','vector'); close(f);
    fprintf(fid,'layers (nm): %s | names: %s\n\n', sprintf('%.1f ',xnm(lsn)), strjoin(string(meta.layerNames),' / '));
  end
end
T = cell2table(rows,'VariableNames',{'Architecture','Condition','V_V','Jterm_uAcm2','JtotSpread_uAcm2','rJ','ElectronResidual','HoleResidual','nMin','pMin','nMax','pMax'});
writetable(T,fullfile(outDir,'spatial_profiles_summary.csv')); disp(T);
end

function [states, J] = continuationSweep(dev, G, Vgrid)
N = dev.N; idx0 = find(abs(Vgrid-0.02)<1e-12,1);
r = gummelSolveQFPAnalytic(dev,0.02,G,struct('psiInit',dev.psiEq,'phinInit',zeros(N,1),'phipInit',zeros(N,1)));
states = cell(numel(Vgrid),1); states{idx0} = r; prev = r;
for v = [0.01,0.005,0]
  r0 = gummelSolveQFPAnalytic(dev,v,G,struct('psiInit',prev.psi,'phinInit',prev.phin,'phipInit',prev.phip));
  prev = r0;
end
states{1} = r0;                                   % V = 0 reached by reverse continuation (grid has no 0.01/0.005 V points)
prev = r;
for k = idx0+1:numel(Vgrid)
  rk = gummelSolveQFPAnalytic(dev,Vgrid(k),G,struct('psiInit',prev.psi,'phinInit',prev.phin,'phipInit',prev.phip));
  states{k} = rk; if rk.isFinite, prev = rk; end
  if rk.isFinite && rk.J < 0 && Vgrid(k) > 0.5, states = states(1:k); break; end
end
J = nan(numel(Vgrid),1);
for k = 1:numel(states), if ~isempty(states{k}) && states{k}.isFinite, J(k) = states{k}.J; end, end
end
