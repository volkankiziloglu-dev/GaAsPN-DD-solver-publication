function generate_heatmap_figures
%GENERATE_HEATMAP_FIGURES Indoor PCE heat map (fig03) and GaP chi-tau heat map (fig04) of the main text.
% Same data, layout and colour maps as generate_main_result_figures.m (sections 'Indoor PCE heatmaps' and
% 'GaP chi-lifetime sensitivity'), but the number in each cell is black or white according to the luminance of
% the cell colour, so it stays readable on light (yellow) and dark cells alike. Only these two figures are
% written; generate_main_result_figures.m also rewrites tables and is not called.
root = fileparts(fileparts(mfilename('fullpath')));
out = fullfile(root,'figures','main_results');
arch = {'HJ-0','F-HJ','RBSF-GaP','RHJ'}; spectra = {'2700K','4000K','6500K','4000F'};
set(groot,'defaultAxesFontName','Arial','defaultAxesFontSize',8,'defaultAxesBox','off','defaultFigureColor','w');

%% Indoor PCE heat maps
I = readtable(fullfile(root,'results','paper_indoor_production','indoor_production_summary_unified.csv'),VariableNamingRule='preserve');
f = figure('Units','inches','Position',[1 1 7.15 4.25]); tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
cmin = min(I.PCE_pct); cmax = max(I.PCE_pct); cm = cividisMap(256);
for il = 1:2
  lx = [100 1000]; ax = nexttile; M = nan(numel(arch),4); C = cell(numel(arch),4);
  for ia = 1:numel(arch), for is = 1:4
    q = I(I.Lux==lx(il) & strcmp(I.Architecture,arch{ia}) & strcmp(I.Spectrum,spectra{is}),:);
    M(ia,is) = q.PCE_pct; C{ia,is} = q.CurveClass{1};
  end, end
  imagesc(ax,M,[cmin cmax]); colormap(ax,cm);
  xticks(ax,1:4); xticklabels(ax,spectra); yticks(ax,1:numel(arch)); yticklabels(ax,arch);
  title(ax,sprintf('(%c) %d lx',char(96+il),lx(il)));
  for ia = 1:numel(arch), for is = 1:4
    suf = ''; if strcmp(C{ia,is},'SCREEN'), suf = '*'; end
    text(ax,is,ia,sprintf('%.1f%s',M(ia,is),suf),'HorizontalAlignment','center','FontSize',7, ...
      'Color',textColor(M(ia,is),cmin,cmax,cm));
  end, end
end
cb = colorbar; cb.Layout.Tile = 'east'; ylabel(cb,'PCE (%)'); saveFig(f,out,'fig03_indoor_pce_comparison');

%% GaP chi-lifetime sensitivity
G = readtable(fullfile(root,'results','paper_gap_sensitivity_am','gap_sensitivity_summary.csv'),VariableNamingRule='preserve');
chis = unique(G.GaPChi_eV)'; taus = unique(G.GaPTau_s)'; gv = G.PCE_pct(~isnan(G.PCE_pct));
lo = min(gv); hi = max(gv); cm = parula(256);
f = figure('Units','inches','Position',[1 1 7.15 5.2]); tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
for ia = 1:numel(arch)
  ax = nexttile; M = nan(numel(taus),numel(chis));
  for it = 1:numel(taus), for ic = 1:numel(chis)
    q = G(strcmp(G.Architecture,arch{ia}) & G.GaPChi_eV==chis(ic) & G.GaPTau_s==taus(it),:);
    if ~isempty(q), M(it,ic) = q.PCE_pct; end
  end, end
  imagesc(ax,M,[lo hi]); set(ax,'YDir','normal'); colormap(ax,cm); hold(ax,'on');
  xticks(ax,1:numel(chis)); xticklabels(ax,string(chis)); yticks(ax,1:numel(taus)); yticklabels(ax,compose('10^{%d}',round(log10(taus))));
  xlabel(ax,'\chi_{GaP} (eV)'); ylabel(ax,'\tau_{GaP} (s)'); title(ax,arch{ia});
  for it = 1:numel(taus), for ic = 1:numel(chis)
    if isnan(M(it,ic))
      plot(ax,[ic-.35 ic+.35],[it-.35 it+.35],'w-','LineWidth',1.2); plot(ax,[ic-.35 ic+.35],[it+.35 it-.35],'w-','LineWidth',1.2);
    else
      text(ax,ic,it,sprintf('%.2f',M(it,ic)),'HorizontalAlignment','center','FontSize',6,'Color',textColor(M(it,ic),lo,hi,cm));
    end
  end, end
end
cb = colorbar(ax); cb.Layout.Tile = 'east'; ylabel(cb,'PCE (%)'); saveFig(f,out,'fig04_gap_chi_lifetime_pce_sensitivity');
end

function c = textColor(v,lo,hi,cm)
% black or white, whichever gives the higher WCAG contrast ratio against the cell colour
k = 1 + round((size(cm,1)-1) * min(max((v-lo)/max(hi-lo,eps),0),1));
rgb = cm(k,:); g = rgb/12.92; h = rgb > 0.04045; g(h) = ((rgb(h)+0.055)/1.055).^2.4;   % sRGB -> linear
L = 0.2126*g(1) + 0.7152*g(2) + 0.0722*g(3);   % relative luminance (WCAG)
if (L+0.05)/0.05 >= 1.05/(L+0.05), c = [0 0 0]; else, c = [1 1 1]; end   % higher contrast ratio wins
end
function saveFig(f,out,name)
exportgraphics(f,fullfile(out,[name '.pdf']),'ContentType','vector'); exportgraphics(f,fullfile(out,[name '.png']),'Resolution',600); close(f)
end
function cm = cividisMap(n)
anchors = [0 .126 .302; .227 .282 .420; .486 .455 .388; .733 .651 .204; 1 .914 .271]; x = linspace(0,1,size(anchors,1)); xi = linspace(0,1,n); cm = interp1(x,anchors,xi);
end
