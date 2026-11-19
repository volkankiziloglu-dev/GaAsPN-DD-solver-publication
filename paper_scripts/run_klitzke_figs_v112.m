function run_klitzke_figs_v112(parts)
%RUN_KLITZKE_FIGS_V112 Klitzke figures and the RHJ doping table with this project's device definition (dd_solver v1.1.2; 2026-11-13).
% Devices as run_klitzke_table_v112.m (AlInP optical table, 85/45 nm MgF2/ZnS ARC; RHJ N_D = 2e16 unless varied, nominal mesh,
% analytic-Jacobian QFP production solver, intrinsic). Replaces the Klitzke_Validasyon scripts plot_klitzke_jv_model_vs_measured.m
% and klitzke_spatial_profiles.m (RHJ N_D = 5e17, other optics) and the doping table of the Supplement.
%   'doping'   RHJ PCE vs absorber doping (2e16 ... 5e17) at 100 and 1000 lx -> rhj_doping_v112.csv
%   'profiles' RHJ at 1000 lx: J_n, J_p, J_tot, n, p at V = 0, MPP and Voc -> fig_klitzke_rhj_profiles.pdf + summary csv
%   'jv'       model (results/paper_klitzke_one_solver, mesh 1) vs digitized measurement -> fig_klitzke_jv_model_vs_measured.pdf
% Output: results/paper_klitzke_figs_v112/ (figures are copied to figures/main_results/ only by hand).
if nargin < 1, parts = {'doping','profiles','jv'}; end
root = fileparts(fileparts(mfilename('fullpath')));
restoredefaultpath; addpath(fullfile(root,'src'),'-begin'); assertCanonicalPath(); warning('off','all');
out = fullfile(root,'results','paper_klitzke_figs_v112'); if ~isfolder(out), mkdir(out); end
alOpt = loadAlInPOptical(); arcFn = @(lam) tmmReflectance(lam,[1,1.38,2.35,cidx(alOpt,lam)],[85,45]);
Tspd = readtable(resolveDataFile('Klitzke2025_CWLED_SPD_digitized.csv')); led = [Tspd.lambda_nm, Tspd.relative_irradiance];
specAt = @(P) [led(:,1), led(:,2)*(P/spectrumPower(led))];

if any(strcmp(parts,'doping'))
  rows = [];
  for Nd = [2e16 5e16 1e17 2e17 5e17]
    dev = makeRHJ(Nd,alOpt);
    for P = [0.31 3.09]
      sp = specAt(P); G = deviceGenAndTransmitML(dev,sp,arcFn);
      [m,cls] = sweep(dev,G,(0:0.02:1.60)',spectrumPower(sp));
      rows(end+1,:) = [Nd, round(P/0.00309), m.Jsc*1e6, m.Voc*1e3, m.FF*100, m.PCE]; %#ok<AGROW>
      fprintf('Nd %.0e %4d lx: Jsc %.3f uA Voc %.1f FF %.2f PCE %.3f %s\n',Nd,rows(end,2),rows(end,3:6),cls);
    end
  end
  writetable(array2table(rows,'VariableNames',{'Nd_cm3','Lux','Jsc_uAcm2','Voc_mV','FF_pct','PCE_pct'}),fullfile(out,'rhj_doping_v112.csv'));
end

if any(strcmp(parts,'profiles'))
  dev = makeRHJ(2e16,alOpt); sp = specAt(3.09); G = deviceGenAndTransmitML(dev,sp,arcFn); N = dev.N;
  ws = @(s) struct('psiInit',s.psi,'phinInit',s.phin,'phipInit',s.phip);
  r02 = gummelSolveQFPAnalytic(dev,0.02,G,struct('psiInit',dev.psiEq,'phinInit',zeros(N,1),'phipInit',zeros(N,1))); prev = r02;
  for v = [0.01 0.005 0], r0 = gummelSolveQFPAnalytic(dev,v,G,ws(prev)); prev = r0; end
  Vg = (0:0.02:1.6)'; st = {r0}; J = r0.J; prev = r0;
  for k = 2:numel(Vg)
    rk = gummelSolveQFPAnalytic(dev,Vg(k),G,ws(prev)); st{end+1} = rk; J(end+1) = rk.J; prev = rk; %#ok<AGROW>
    if rk.J < 0, break; end
  end
  [~,km] = max(Vg(1:numel(J)).*J(:).*(J(:)>0));
  kpos = find(J > 0,1,'last'); V1 = Vg(kpos); J1 = J(kpos); V2 = Vg(kpos+1); J2 = J(kpos+1); prev = st{kpos};
  for it = 1:8
    Vo = V1 - J1*(V2-V1)/(J2-J1); rOc = gummelSolveQFPAnalytic(dev,Vo,G,ws(prev));
    if abs(rOc.J) < 1e-4*abs(J(1)), break; end
    if rOc.J > 0, V1 = Vo; J1 = rOc.J; prev = rOc; else, V2 = Vo; J2 = rOc.J; end
  end
  sel = {st{1}, st{km}, rOc}; Vsel = [0, Vg(km), Vo];
  labs = {'V = 0 (J_{sc})', sprintf('V = %.2f V (MPP)',Vg(km)), sprintf('V = %.3f V (V_{oc})',Vo)};
  xnm = dev.x*1e7; xseg = 0.5*(xnm(1:end-1)+xnm(2:end)); lsn = dev.layerStartNode; rows = {};
  f = figure('Visible','off','Units','centimeters','Position',[2 2 24 14],'Color','w');
  tl = tiledlayout(f,2,3,'TileSpacing','compact','Padding','compact');
  for k = 1:3
    s = sel{k}; Jn = s.JnSeg*1e6; Jp = s.JpSeg*1e6; Jt = Jn + Jp; sp_ = max(Jt) - min(Jt);
    ax = nexttile(tl,k); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
    for j = 2:numel(lsn), xline(ax,xnm(lsn(j))/1e3,':','Color',[.55 .55 .55],'HandleVisibility','off'); end
    yline(ax,0,'-','Color',[.8 .8 .8],'HandleVisibility','off');
    h1 = plot(ax,xseg/1e3,Jn,'-','Color',[0 .447 .741],'LineWidth',1.4); h2 = plot(ax,xseg/1e3,Jp,'-','Color',[.85 .325 .098],'LineWidth',1.4);
    h3 = plot(ax,xseg/1e3,Jt,'k--','LineWidth',1.1); title(ax,labs{k},'FontWeight','normal'); xlim(ax,[0 xnm(end)/1e3]);
    if k == 1, ylabel(ax,'J (\muA cm^{-2})'); legend(ax,[h1 h2 h3],{'J_n','J_p','J_{tot}'},'Location','east','Box','off'); end
    ax = nexttile(tl,k+3); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
    for j = 2:numel(lsn), xline(ax,xnm(lsn(j))/1e3,':','Color',[.55 .55 .55],'HandleVisibility','off'); end
    g1 = plot(ax,xnm/1e3,s.n,'-','Color',[0 .447 .741],'LineWidth',1.4); g2 = plot(ax,xnm/1e3,s.p,'-','Color',[.85 .325 .098],'LineWidth',1.4);
    set(ax,'YScale','log'); xlim(ax,[0 xnm(end)/1e3]); ylim(ax,[1e0 1e20]); yticks(ax,10.^(0:4:20)); xlabel(ax,'x (\mum)');
    if k == 1, ylabel(ax,'n, p (cm^{-3})'); legend(ax,[g1 g2],{'n','p'},'Location','southeast','Box','off'); end
    rows(end+1,:) = {Vsel(k), s.J*1e6, sp_, sp_/max(abs(s.J*1e6),1e-9), s.currentConservation.rJ, ...
      s.continuityResiduals.electron.maxRelGlobal, s.continuityResiduals.hole.maxRelGlobal}; %#ok<AGROW>
  end
  title(tl,'Klitzke lower-doped RHJ (N_D = 2\times10^{16} cm^{-3}), 1000 lx (CW-LED), intrinsic model','FontWeight','normal');
  exportgraphics(f,fullfile(out,'fig_klitzke_rhj_profiles.png'),'Resolution',300); exportgraphics(f,fullfile(out,'fig_klitzke_rhj_profiles.pdf'),'ContentType','vector'); close(f);
  T = cell2table(rows,'VariableNames',{'V_V','Jterm_uAcm2','JtotSpread_uAcm2','JtotSpreadRel','solver_rJ','ElectronResidual','HoleResidual'});
  writetable(T,fullfile(out,'klitzke_rhj_1000lx_profiles_summary.csv')); disp(T);
end

if any(strcmp(parts,'jv'))
  Tm = readtable(resolveDataFile(fullfile('experimental','klitzke2025','klitzke2025_fig4_iv_curves.csv')));
  VocT = struct('HJ',[1.04 1.24],'RHJ',[1.28 1.35]);   % Klitzke Table II, 100 / 1000 lx
  lxs = [1000 100]; devs = {'Homojunction','LowerRHJ'}; names = {'HJ','RHJ'};
  cHJ = [0 119 187]/255; cRHJ = [204 51 17]/255; cols = {cHJ,cRHJ};
  f = figure('Visible','off','Units','inches','Position',[1 1 7.15 3.6],'Color','w');
  tiledlayout(1,2,'TileSpacing','compact','Padding','compact'); top = 0;
  for k = 1:2
    lx = lxs(k); ax = nexttile; hold(ax,'on'); box(ax,'on'); grid(ax,'on');
    for d = 1:2
      nm = names{d}; c = cols{d};
      m = readtable(fullfile(root,'results','paper_klitzke_one_solver',sprintf('%s_%dlx_mesh1_qfp.csv',nm,lx)));
      m = m(m.Finite==1,:); mv = m.V_V; mj = m.J_Acm2*1e6; ok = mj > -0.05*max(mj);
      plot(ax,mv(ok),mj(ok),'-','Color',c,'LineWidth',1.5,'DisplayName',[nm ' model']); top = max(top,max(mj));
      q = Tm(strcmp(Tm.device,devs{d}) & Tm.lux==lx,:); [vq,j] = sort(q.V_V); jq = q.J_uAcm2(j);
      ks = max(1,round(numel(vq)/70)); plot(ax,vq(1:ks:end),jq(1:ks:end),'o','Color',c,'MarkerSize',3.2,'MarkerFaceColor','w','DisplayName',[nm ' measured']);
      xline(ax,VocT.(nm)(1+(lx==1000)),':','Color',c,'HandleVisibility','off');
    end
    yline(ax,0,'k-','HandleVisibility','off'); xlim(ax,[0 1.6]); ylim(ax,[-0.04 1.12]*top); top = 0;
    xlabel(ax,'Voltage (V)'); if k==1, ylabel(ax,'Current density (\muA cm^{-2})'); end
    title(ax,sprintf('(%c) %d lx',char(96+k),lx),'HorizontalAlignment','left');
    if k==1, legend(ax,'Location','southwest','Box','off','FontSize',7); end
  end
  exportgraphics(f,fullfile(out,'fig_klitzke_jv_model_vs_measured.png'),'Resolution',300);
  exportgraphics(f,fullfile(out,'fig_klitzke_jv_model_vs_measured.pdf'),'ContentType','vector'); close(f);
end
end

function [m,cls] = sweep(dev,G,Vg,Pin)
N = dev.N; i2 = find(abs(Vg-0.02)<1e-12,1); s = @gummelSolveQFPAnalytic;
r02 = s(dev,0.02,G,struct('psiInit',dev.psiEq,'phinInit',zeros(N,1),'phipInit',zeros(N,1))); prev = r02;
for v = [0.01 0.005 0], r0 = s(dev,v,G,struct('psiInit',prev.psi,'phinInit',prev.phin,'phipInit',prev.phip)); prev = r0; end
S = repmat(r02,numel(Vg),1); S(1) = r0; S(i2) = r02; prev = r02;
for k = i2+1:numel(Vg), rk = s(dev,Vg(k),G,struct('psiInit',prev.psi,'phinInit',prev.phin,'phipInit',prev.phip)); S(k) = rk; if rk.isFinite, prev = rk; end, end
g = indoorAcceptanceGate(dev,G,Vg,S,Pin); f = g.isFinite(:); m = pvMetrics([g.Vs(f), g.Jterm(f)],Pin); cls = char(g.curveClass);
end

function dev = makeRHJ(Nd,alOpt)
L(1)=layer(25e-7,30,2.35,3.78,11.8,8e18,0,100,10,1e-9,2.5e18,7e18,0,.05); L(1).useOpticalTable=true; L(1).opticalTable=alOpt;
L(2)=layer(850e-7,190,1.90,4.10,11.8,Nd,0,500,30,248e-9,6.5e17,1.5e19,3e4,.05); L(2).Brad=1e-10; L(2).srhCenters=struct('tauN',337e-9,'tauP',248e-9,'EtRelEi_eV',0);
L(3)=layer(100e-7,30,2.35,3.78,11.8,0,2e17,100,10,1e-9,2.5e18,7e18,0,.05);
dev = setupDeviceML(L,2,struct('Sn_int',{0,0},'Sp_int',{0,0}),struct('Sback',1e5));
end
function z=cidx(opt,lam)
E=1239.842/lam; n=interp1(opt.E,opt.n,E,'linear','extrap'); a=interp1(opt.E,opt.alpha,E,'linear',0); k=max(a,0)*(lam*1e-7)/(4*pi); z=n-1i*k;
end
function L=layer(t,N,Eg,chi,eps_,Nd,Na,muN,muP,tau,Nc,Nv,alphaMax,edgeWidth)
L=struct('t',t,'N',N,'Eg',Eg,'chi',chi,'eps',eps_,'Nd',Nd,'Na',Na,'muN',muN,'muP',muP,'tau',tau,'Nc',Nc,'Nv',Nv,'NcNv',sqrt(Nc*Nv),'Brad',0,'CnAuger',3e-30,'CpAuger',3e-30,'alphaMax',alphaMax,'edgeWidth',edgeWidth,'useOpticalTable',false,'opticalTable',[],'srhCenters',[]);
end
