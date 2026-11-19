% Single-column version of fig01_outdoor_jv_comparison (AM1.5G J-V, SRH+Rad+Auger) for manuscript_en_kisa.
% Same data and styling as generate_main_result_figures.m; output fig01_outdoor_jv_comparison_1col.{pdf,png}.
root = fileparts(fileparts(mfilename('fullpath')));
out = fullfile(root,'figures','main_results');
arch = {'HJ-0','F-HJ','RBSF-GaP','RHJ'}; names = {'HJ-0','F-HJ','RBSF–GaP','RHJ'};
cols = [0 119 187;238 119 51;204 51 17;102 102 102]/255; marks = {'o','d','v','x'};
f = figure('Units','inches','Position',[1 1 3.5 2.7],'Color','w'); ax = axes(f); hold(ax,'on');
for i = 1:numel(arch)
  D = readtable(fullfile(root,'results','paper_recombination_scenarios',sprintf('SRH+Rad+Auger_%s_jv.csv',arch{i})));
  ok = isfinite(D.J_Acm2);
  plot(ax,D.V_V(ok),1e3*D.J_Acm2(ok),'Color',cols(i,:),'LineWidth',1.1,'Marker',marks{i}, ...
    'MarkerIndices',unique(round(linspace(1,nnz(ok),10))),'MarkerSize',3,'DisplayName',names{i});
end
yline(ax,0,'k-','HandleVisibility','off'); xlim(ax,[0 1.58]); ylim(ax,[-1.5 18]); grid(ax,'on');
set(ax,'FontName','Arial','FontSize',7,'LineWidth',0.6);
xlabel(ax,'Voltage (V)'); ylabel(ax,'Current density (mA cm^{-2})');
legend(ax,'Location','west','NumColumns',2,'Box','off','FontSize',7);
exportgraphics(f,fullfile(out,'fig01_outdoor_jv_comparison_1col.pdf'),'ContentType','vector');
exportgraphics(f,fullfile(out,'fig01_outdoor_jv_comparison_1col.png'),'Resolution',600); close(f);
