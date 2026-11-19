% PCE versus illuminance (CIE 6500 K, 50-5000 lx) of the four architectures, intrinsic curves
% (results/paper_voc_lux/lux_series_pce.csv from cases/analyze_lux_series.py). Filled markers: PASS;
% open markers: SCREEN class of the residual diagnostic. Output: figures/main_results/fig_lux_pce.{pdf,png}
root = fileparts(fileparts(mfilename('fullpath')));
T = readtable(fullfile(root,'results','paper_voc_lux','lux_series_pce.csv'));
arch = {'HJ-0','F-HJ','RBSF-GaP','RHJ'}; names = {'HJ-0','F-HJ','RBSF–GaP','RHJ'};
cols = [0 119 187;238 119 51;204 51 17;102 102 102]/255; marks = {'o','d','v','s'};
f = figure('Units','inches','Position',[1 1 3.5 2.7],'Color','w'); ax = axes(f); hold(ax,'on'); box(ax,'on');
for i = 1:4
  S = T(strcmp(T.Architecture,arch{i}),:); S = sortrows(S,'lux');
  plot(ax,S.lux,S.PCE_pct,'-','Color',cols(i,:),'LineWidth',1.1,'DisplayName',names{i});
  scr = strcmp(S.Class,'SCREEN');
  plot(ax,S.lux(~scr),S.PCE_pct(~scr),marks{i},'Color',cols(i,:),'MarkerFaceColor',cols(i,:),'MarkerSize',3.5,'HandleVisibility','off');
  plot(ax,S.lux(scr),S.PCE_pct(scr),marks{i},'Color',cols(i,:),'MarkerFaceColor','w','MarkerSize',3.5,'HandleVisibility','off');
end
set(ax,'XScale','log','FontName','Arial','FontSize',7,'LineWidth',0.6); grid(ax,'on');
xlim(ax,[40 6000]); ylim(ax,[36 49]); xticks(ax,[50 100 200 500 1000 2000 5000]);
xticklabels(ax,{'50','100','200','500','1000','2000','5000'});
xlabel(ax,'Illuminance (lx)'); ylabel(ax,'PCE (%)');
legend(ax,'Location','northwest','Box','off','FontSize',7);
exportgraphics(f,fullfile(root,'figures','main_results','fig_lux_pce.pdf'),'ContentType','vector');
exportgraphics(f,fullfile(root,'figures','main_results','fig_lux_pce.png'),'Resolution',600); close(f);
