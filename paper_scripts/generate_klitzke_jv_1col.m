function generate_klitzke_jv_1col(variant)
% Single-column versions of fig_klitzke_jv_model_vs_measured (Fig. 4 of manuscript_en_kisa).
% Same data and styling as run_klitzke_figs_v112('jv').
%   'stacked'    (a) 1000 lx above (b) 100 lx, shared voltage axis -> fig_klitzke_jv_1col_stacked.{pdf,png}
%   'normalized' one panel, J/Jsc(model) at both illuminances    -> fig_klitzke_jv_1col_norm.{pdf,png}
if nargin < 1, variant = 'stacked'; end
root = fileparts(fileparts(mfilename('fullpath')));
out = fullfile(root,'figures','main_results');
Tm = readtable(resolveDataFile(fullfile('experimental','klitzke2025','klitzke2025_fig4_iv_curves.csv')));
VocT = struct('HJ',[1.04 1.24],'RHJ',[1.28 1.35]);   % Klitzke Table II, 100 / 1000 lx
lxs = [1000 100]; devs = {'Homojunction','LowerRHJ'}; names = {'HJ','RHJ'};
cols = {[0 119 187]/255,[204 51 17]/255};
switch variant
  case 'stacked'
    f = figure('Visible','off','Units','inches','Position',[1 1 3.5 4.6],'Color','w');
    tl = tiledlayout(f,2,1,'TileSpacing','tight','Padding','compact');
    for k = 1:2
      lx = lxs(k); ax = nexttile(tl); hold(ax,'on'); box(ax,'on'); grid(ax,'on'); top = 0;
      for d = 1:2
        [mv,mj,vq,jq] = getData(root,Tm,names{d},devs{d},lx); c = cols{d}; top = max(top,max(mj));
        plot(ax,mv,mj,'-','Color',c,'LineWidth',1.2,'DisplayName',[names{d} ' model']);
        plot(ax,vq,jq,'o','Color',c,'MarkerSize',2.6,'MarkerFaceColor','w','DisplayName',[names{d} ' measured']);
        xline(ax,VocT.(names{d})(1+(lx==1000)),':','Color',c,'HandleVisibility','off');
      end
      yline(ax,0,'k-','HandleVisibility','off'); xlim(ax,[0 1.5]); ylim(ax,[-0.04 1.15]*top);
      set(ax,'FontName','Arial','FontSize',7,'LineWidth',0.6);
      ylabel(ax,'J (\muA cm^{-2})');
      text(ax,0.03,0.08,sprintf('(%c) %d lx',char(96+k),lx),'Units','normalized','FontSize',8,'FontName','Arial');
      if k == 1, xticklabels(ax,{}); legend(ax,'Location','southwest','Position',[0.2 0.6 0.3 0.1],'Box','off','FontSize',6.5,'NumColumns',1);
      else, xlabel(ax,'Voltage (V)'); end
    end
    base = 'fig_klitzke_jv_1col_stacked';
  case 'overlay'
    % one panel, absolute J on a log axis so that 100 and 1000 lx share the axis
    f = figure('Visible','off','Units','inches','Position',[1 1 3.5 3.0],'Color','w');
    ax = axes(f); hold(ax,'on'); box(ax,'on'); grid(ax,'on'); ls = {'-','--'}; mk = {'o','s'};
    for k = 1:2
      for d = 1:2
        [mv,mj,vq,jq] = getData(root,Tm,names{d},devs{d},lxs(k)); c = cols{d};
        p = mj > 0; plot(ax,mv(p),mj(p),ls{k},'Color',c,'LineWidth',1.1,'DisplayName',sprintf('%s model',names{d}));
        p = jq > 0; plot(ax,vq(p),jq(p),mk{k},'Color',c,'MarkerSize',2.4,'MarkerFaceColor','w','DisplayName',sprintf('%s measured',names{d}));
        xline(ax,VocT.(names{d})(1+(lxs(k)==1000)),':','Color',c,'HandleVisibility','off');
      end
    end
    set(ax,'YScale','log','FontName','Arial','FontSize',7,'LineWidth',0.6,'YMinorGrid','off','YMinorTick','off');
    xlim(ax,[0 1.45]); ylim(ax,[1 200]); yticks(ax,[1 2 5 10 20 50 100 200]);
    xlabel(ax,'Voltage (V)'); ylabel(ax,'J (\muA cm^{-2})');
    text(ax,0.05,140,'1000 lx (solid, circles)','FontSize',6.5,'FontName','Arial');
    text(ax,0.05,14,'100 lx (dashed, squares)','FontSize',6.5,'FontName','Arial');
    h = findobj(ax,'Type','line'); h = flipud(h); h = h([1 2 3 4]);
    legend(ax,h,'Location','southwest','Box','off','FontSize',6.5);
    base = 'fig_klitzke_jv_1col_overlay';
  case 'overlaylin'
    % one panel, absolute J on a linear axis
    f = figure('Visible','off','Units','inches','Position',[1 1 3.5 3.0],'Color','w');
    ax = axes(f); hold(ax,'on'); box(ax,'on'); grid(ax,'on'); ls = {'-','--'}; mk = {'o','s'};
    for k = 1:2
      for d = 1:2
        [mv,mj,vq,jq] = getData(root,Tm,names{d},devs{d},lxs(k)); c = cols{d};
        plot(ax,mv,mj,ls{k},'Color',c,'LineWidth',1.1,'DisplayName',sprintf('%s model',names{d}));
        plot(ax,vq,jq,mk{k},'Color',c,'MarkerSize',2.4,'MarkerFaceColor','w','DisplayName',sprintf('%s measured',names{d}));
      end
    end
    yline(ax,0,'k-','HandleVisibility','off');
    set(ax,'FontName','Arial','FontSize',7,'LineWidth',0.6); xlim(ax,[0 1.45]); ylim(ax,[-3 125]);
    xlabel(ax,'Voltage (V)'); ylabel(ax,'J (\muA cm^{-2})');
    text(ax,0.05,119,'1000 lx (solid, circles)','FontSize',6.5,'FontName','Arial');
    text(ax,0.05,17,'100 lx (dashed, squares)','FontSize',6.5,'FontName','Arial');
    h = flipud(findobj(ax,'Type','line')); h = h(1:4);
    legend(ax,h,'Location','west','Box','off','FontSize',6.5);
    base = 'fig_klitzke_jv_1col_overlay_lin';
  case 'normalized'
    f = figure('Visible','off','Units','inches','Position',[1 1 3.5 2.8],'Color','w');
    ax = axes(f); hold(ax,'on'); box(ax,'on'); grid(ax,'on'); ls = {'-','--'}; mk = {'o','s'};
    for k = 1:2
      for d = 1:2
        [mv,mj,vq,jq] = getData(root,Tm,names{d},devs{d},lxs(k)); c = cols{d}; js = mj(1);
        plot(ax,mv,mj/js,ls{k},'Color',c,'LineWidth',1.1,'DisplayName',sprintf('%s model, %d lx',names{d},lxs(k)));
        plot(ax,vq,jq/js,mk{k},'Color',c,'MarkerSize',2.4,'MarkerFaceColor','w','DisplayName',sprintf('%s measured, %d lx',names{d},lxs(k)));
      end
    end
    yline(ax,0,'k-','HandleVisibility','off'); xlim(ax,[0 1.5]); ylim(ax,[-0.03 1.1]);
    set(ax,'FontName','Arial','FontSize',7,'LineWidth',0.6);
    xlabel(ax,'Voltage (V)'); ylabel(ax,'J / J_{sc,model}');
    legend(ax,'Location','southwest','Box','off','FontSize',5.5,'NumColumns',1);
    base = 'fig_klitzke_jv_1col_norm';
end
exportgraphics(f,fullfile(out,[base '.pdf']),'ContentType','vector');
exportgraphics(f,fullfile(out,[base '.png']),'Resolution',300); close(f);
end

function [mv,mj,vq,jq] = getData(root,Tm,nm,dev,lx)
if strcmp(nm,'HJ')   % homojunction with the back-surface field of Klitzke et al. Fig. 1(a) (run_klitzke_hj_bsf.m)
  m = readtable(fullfile(root,'results','paper_klitzke_hj_bsf',sprintf('HJBSF_AlInP_%d_curve.csv',lx)));
else
  m = readtable(fullfile(root,'results','paper_klitzke_one_solver',sprintf('%s_%dlx_mesh1_qfp.csv',nm,lx)));
end
m = m(m.Finite==1,:); mv = m.V_V; mj = m.J_Acm2*1e6; ok = mj > -0.05*max(mj); mv = mv(ok); mj = mj(ok);
q = Tm(strcmp(Tm.device,dev) & Tm.lux==lx,:); [vq,j] = sort(q.V_V); jq = q.J_uAcm2(j);
ks = max(1,round(numel(vq)/70)); vq = vq(1:ks:end); jq = jq(1:ks:end);
end
