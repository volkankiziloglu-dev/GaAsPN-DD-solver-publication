function generate_spatial_comparison_figure()
%GENERATE_SPATIAL_COMPARISON_FIGURE HJ-0 vs RBSF-GaP at the AM1.5G maximum power point: minority-electron density n(x) and electron current J_n(x)
% from the profiles of run_paper_spatial_profiles.m (results/paper_spatial_profiles). Output: figures/main_results/fig_spatial_comparison_MPP.{png,pdf}.
root = fileparts(fileparts(mfilename('fullpath'))); d = fullfile(root,'results','paper_spatial_profiles');
A = {'HJ-0','RBSF-GaP'}; col = {[.45 .45 .45],[0 .447 .741]}; ls = {'-','-'};
f = figure('Visible','off','Units','centimeters','Position',[2 2 17 7],'Color','w'); tl = tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');
ax1 = nexttile(tl); hold(ax1,'on'); box(ax1,'on'); grid(ax1,'on'); ax2 = nexttile(tl); hold(ax2,'on'); box(ax2,'on'); grid(ax2,'on');
for k = 1:2
  c = readtable(fullfile(d,sprintf('%s_AM1.5G_MPP_carriers.csv',A{k}))); j = readtable(fullfile(d,sprintf('%s_AM1.5G_MPP_currents.csv',A{k})));
  plot(ax1,c.x_nm/1e3,c.n_cm3,ls{k},'Color',col{k},'LineWidth',1.5,'DisplayName',A{k});
  plot(ax2,j.x_nm/1e3,j.Jn_uAcm2/1e3,ls{k},'Color',col{k},'LineWidth',1.5,'DisplayName',A{k});
end
set(ax1,'YScale','log'); ylim(ax1,[1e8 1e20]); xlim(ax1,[0 1.62]); xlim(ax2,[0 1.62]);
xlabel(ax1,'x (\mum)'); ylabel(ax1,'Electron density n (cm^{-3})'); xlabel(ax2,'x (\mum)'); ylabel(ax2,'J_n (mA cm^{-2})');
title(ax1,'(a) minority electrons, MPP','FontWeight','normal'); title(ax2,'(b) electron current, MPP','FontWeight','normal');
legend(ax1,'Location','southwest','Box','off');
out = fullfile(root,'figures','main_results'); exportgraphics(f,fullfile(out,'fig_spatial_comparison_MPP.png'),'Resolution',300); exportgraphics(f,fullfile(out,'fig_spatial_comparison_MPP.pdf'),'ContentType','vector'); close(f);
end
