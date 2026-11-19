function plotPaperStageAJV(runs,pngFile,pdfFile)
%PLOTPAPERSTAGEAJV Publication-scale J-V plot for Stage A results.
%
% Deep-forward diode current can be orders of magnitude larger than the
% photovoltaic current. The displayed limits are therefore tied to Jsc and
% Voc so that the power-producing region remains legible; full numerical
% curves remain available in the exported CSV files.

fig=figure('Visible','off','Color','w','Position',[100 100 1000 680]);
ax=axes(fig); hold(ax,'on');
colors=lines(numel(runs));
styles={'-','--','-.',':','-','--','-.'};
jscMax=0; vocMax=0;
for k=1:numel(runs)
  c=runs(k).curve;
  ok=c.Accepted & isfinite(c.J_Acm2);
  plot(ax,c.V_V(ok),c.J_Acm2(ok)*1e3,'LineWidth',1.8, ...
    'Color',colors(k,:),'LineStyle',styles{k},'DisplayName',runs(k).code);
  jscMax=max(jscMax,runs(k).metrics.Jsc*1e3);
  vocMax=max(vocMax,runs(k).metrics.Voc);
end
yline(ax,0,'k:','HandleVisibility','off');
xlabel(ax,'Voltage (V)'); ylabel(ax,'Current density (mA cm^{-2})');
title(ax,'Controlled AM1.5G comparison of GaAsPN architectures');
legend(ax,'Location','southwest','NumColumns',2);
grid(ax,'on'); box(ax,'on');
xlim(ax,[0,min(1.70,vocMax+0.08)]);
ylim(ax,[-0.15*jscMax,1.10*jscMax]);
set(ax,'FontName','Arial','FontSize',11,'LineWidth',0.8);
exportgraphics(fig,pngFile,'Resolution',300);
exportgraphics(fig,pdfFile,'ContentType','vector');
close(fig);
end
