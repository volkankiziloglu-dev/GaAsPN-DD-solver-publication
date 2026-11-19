% Indoor spectra at 1000 lx (CIE LED-B1/B3/B5, FL11) and the Klitzke CW-LED
% (scaled to 3.09 W/m2), with the photopic V(lambda) and absorber edges.
% Panel (a): LEDs + CW-LED; panel (b): FL11. Output:
% figures/main_results/fig_indoor_spectra.{pdf,png}
here = fileparts(mfilename('fullpath')); root = fileparts(here);
addpath(fullfile(root,'src'));
out = fullfile(root,'figures','main_results');
V = readmatrix(resolveDataFile('CIE_sle_photopic.csv'));
edges = [1239.84/1.8565, 1239.84/1.90, 1239.84/1.9211];   % GaAsPN 3% N, GaInP, GaAsPN 2.5% N
f = figure('Units','inches','Position',[1 1 3.5 3.6],'Color','w');
t = tiledlayout(f, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
Pin = zeros(1,5);
% ---- (a) LEDs ----
ax = nexttile(t); hold(ax,'on'); box(ax,'on');
modes = ["2700K","4000K","6500K"]; names = ["LED-B1 (2700 K)","LED-B3 (4000 K)","LED-B5 (6500 K)"];
cols = [0.85 0.33 0.10; 0.47 0.67 0.19; 0 0.45 0.74];
for k = 1:3
  pts = buildCIEIndoorSpectrum(modes(k), 1000); Pin(k) = trapz(pts(:,1), pts(:,2));
  plot(ax, pts(:,1), 1e3*pts(:,2), 'LineWidth', 1.0, 'Color', cols(k,:), 'DisplayName', names(k));
end
kl = readmatrix(resolveDataFile('Klitzke2025_CWLED_SPD_digitized.csv'));
wl = kl(:,1); s = max(kl(:,2),0); s = s / trapz(wl, s) * 3.09; Pin(5) = trapz(wl, s);
plot(ax, wl, 1e3*s, '-.', 'LineWidth', 1.0, 'Color', [0.3 0.3 0.3], 'DisplayName', 'CW-LED (Klitzke)');
decorate(ax, V, edges, 40, true);
legend(ax, 'Location', 'northeast', 'FontSize', 6, 'Box', 'off');
text(ax, 0.02, 0.93, '(a)', 'Units', 'normalized', 'FontSize', 8);
% ---- (b) FL11 ----
ax = nexttile(t); hold(ax,'on'); box(ax,'on');
pts = buildCIEIndoorSpectrum("4000F", 1000); Pin(4) = trapz(pts(:,1), pts(:,2));
plot(ax, pts(:,1), 1e3*pts(:,2), 'LineWidth', 0.9, 'Color', [0.49 0.18 0.56], 'DisplayName', 'FL11 (4000F)');
decorate(ax, V, edges, 200, false);
legend(ax, 'Location', 'northeast', 'FontSize', 6, 'Box', 'off');
text(ax, 0.02, 0.93, '(b)', 'Units', 'normalized', 'FontSize', 8);
xlabel(t, 'Wavelength (nm)', 'FontSize', 7);
ylabel(t, 'I_\lambda at 1000 lx (mW m^{-2} nm^{-1})', 'FontSize', 7);
exportgraphics(f, fullfile(out,'fig_indoor_spectra.pdf'), 'ContentType', 'vector');
exportgraphics(f, fullfile(out,'fig_indoor_spectra.png'), 'Resolution', 600); close(f);
fprintf('Pin (W/m2) at 1000 lx [B1 B3 B5 FL11 CW]: %s\n', mat2str(Pin, 4));

function decorate(ax, V, edges, ymax, labels)
plot(ax, V(:,1), 0.95*ymax*V(:,2), ':', 'LineWidth', 1.0, 'Color', [0.55 0.55 0.55], ...
     'DisplayName', 'V(\lambda), scaled');
% band edges as short segments (identified by line style in the caption): 3% N dashed, GaInP grey dashed, 2.5% N dash-dotted
h = 0.48*ymax;
plot(ax, edges(1)*[1 1], [0 h], '--', 'Color', [0.2 0.2 0.2], 'LineWidth', 0.8, 'HandleVisibility', 'off');
plot(ax, edges(2)*[1 1], [0 h], '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.8, 'HandleVisibility', 'off');
plot(ax, edges(3)*[1 1], [0 h], '-.', 'Color', [0.2 0.2 0.2], 'LineWidth', 0.8, 'HandleVisibility', 'off');
xlim(ax, [380 780]); ylim(ax, [0 ymax]);
set(ax, 'FontSize', 7, 'LineWidth', 0.6);
end
