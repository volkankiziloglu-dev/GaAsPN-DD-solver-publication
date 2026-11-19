function T=generate_rear_band_alignment_figure()
%GENERATE_REAR_BAND_ALIGNMENT_FIGURE Full-device equilibrium band diagrams.
assertCanonicalPath();root=fileparts(fileparts(mfilename('fullpath')));
out=fullfile(root,'figures','main_results');if ~isfolder(out),mkdir(out);end
codes=["HJ-0","F-HJ","RBSF-GaP","RHJ"];
titles=["HJ-0 reference","Front heterojunction","GaP rear BSF","Rear heterojunction"];
f=figure('Units','inches','Position',[1 1 7.15 5.2],'Color','w');
tl=tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
rows=cell(0,6);
for k=1:numel(codes)
 [dev,meta]=buildPaperArchitecture(codes(k),solveEquilibrium=true);
 xnm=(dev.x-dev.x(1))*1e7; Ec=-dev.psiEq-(dev.chi-dev.chiRef); Ev=Ec-dev.Eg;
 ax=nexttile;hold(ax,'on');
 plot(ax,xnm,Ec,'LineWidth',1.45,'Color',[0 0.447 0.741]);
 plot(ax,xnm,Ev,'LineWidth',1.45,'Color',[0.85 0.325 0.098]);
 yline(ax,0,':','E_F','Color',[.25 .25 .25],'LabelHorizontalAlignment','left','FontSize',7);
 if numel(meta.layers)>1
  for j=2:numel(meta.layers)
   left=meta.layers(j-1); right=meta.layers(j);
   isPhysical = string(left.material)~=string(right.material) || ...
       left.Nd~=right.Nd || left.Na~=right.Na;
   if isPhysical
    xj=sum([meta.layers(1:j-1).t])*1e7;
    xline(ax,xj,'--','Color',[.55 .55 .55],'LineWidth',.6);
   end
  end
 end
 grid(ax,'on');box(ax,'on');xlim(ax,[0 1700]);ylim(ax,[-2.2 2.5]);title(ax,titles(k),'FontWeight','normal');
 if k>=3,xlabel(ax,'Position from front surface (nm)');end
 if mod(k,2)==1,ylabel(ax,'Energy (eV)');end
 if k==1,legend(ax,{'E_C','E_V','E_F'},'Location','southwest','Box','off');end
 rear=meta.layers(end);prev=meta.layers(end-1);
 rows(end+1,:)={codes(k),prev.material,rear.material,rear.Na,rear.chi-prev.chi,(rear.chi+rear.Eg)-(prev.chi+prev.Eg)}; %#ok<AGROW>
end
title(tl,'Full-device equilibrium band alignment at 300 K and V_{app}=0','FontWeight','normal');
exportgraphics(f,fullfile(out,'fig07_rear_band_alignment.pdf'),'ContentType','vector');
exportgraphics(f,fullfile(out,'fig07_rear_band_alignment.png'),'Resolution',600);close(f);
T=cell2table(rows,'VariableNames',{'Architecture','AbsorberMaterial','RearMaterial','RearNa_cm3','DeltaChi_eV','DeltaChiPlusEg_eV'});
writetable(T,fullfile(out,'table_rear_band_offsets.csv'));disp(T);
end
