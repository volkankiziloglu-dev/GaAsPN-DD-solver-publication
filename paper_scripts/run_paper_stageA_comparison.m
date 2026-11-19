function summary = run_paper_stageA_comparison()
%RUN_PAPER_STAGEA_COMPARISON Controlled AM1.5G comparison for the paper.
%
% Builds all seven devices through buildPaperArchitecture, applies one common
% MgF2/ZnS ARC, solves one common voltage grid, filters every metric through
% qfpResultAccepted, and exports the comparison table and J-V figure.

assertCanonicalPath();
warning('off','MATLAB:nearlySingularMatrix');
warning('off','MATLAB:singularMatrix');
root=fileparts(fileparts(mfilename('fullpath')));
resultsDir=fullfile(root,'results','paper_stageA');
reportsDir=fullfile(root,'reports');
if ~exist(resultsDir,'dir'), mkdir(resultsDir); end
if ~exist(reportsDir,'dir'), mkdir(reportsDir); end

codes=["HJ-0","HJ-U","HJ-G","F-HJ","RBSF-GaP","RBSF-GaNP","RHJ"];
am15=loadRealAM15G();
Pin=spectrumPower(am15);
[d1nm,d2nm,Rbar]=optimizeCommonARC(am15);
arcFn=@(lam) tmmReflectance(lam, ...
  [1.0,1.38,2.35,gapRefractiveIndexReal(lam)],[d1nm,d2nm]);

fprintf('Stage A common ARC: MgF2=%.1f nm, ZnS=%.1f nm, weighted R=%.5f\n', ...
  d1nm,d2nm,Rbar);

Vgrid=(0:0.02:1.70)';
runs=struct('code',{},'meta',{},'curve',{},'metrics',{},'diagnostics',{});
rows=cell(numel(codes),13);

for a=1:numel(codes)
  code=codes(a);
  fprintf('\n=== %s ===\n',code);
  [dev,meta]=buildPaperArchitecture(code);
  G=deviceGenAndTransmitML(dev,am15,arcFn);
  [curve,diag]=strictSweep(dev,G,Vgrid);
  acceptedCurve=curve(curve.Accepted,["V_V","J_Acm2"]);

  if height(acceptedCurve)>=4 && any(abs(acceptedCurve.V_V)<1e-12) && ...
      any(acceptedCurve.J_Acm2<0)
    m=pvMetrics(table2array(acceptedCurve),Pin);
  else
    m=struct('Jsc',NaN,'Voc',NaN,'FF',NaN,'PCE',NaN,'Pmax',NaN);
  end

  runs(a)=struct('code',code,'meta',meta,'curve',curve,'metrics',m, ...
    'diagnostics',diag); %#ok<AGROW>
  rows(a,:)={code,meta.architectureName,m.Jsc*1e3,m.Voc*1e3,m.FF*100,m.PCE, ...
    diag.nAccepted,diag.nAttempted,diag.maxDeltaJabs,diag.maxDeltaJrel, ...
    diag.maxElectronResidual,diag.maxHoleResidual,diag.maxPoissonResidual};

  writetable(curve,fullfile(resultsDir,code+"_jv.csv"));
  fprintf('%s: accepted %d/%d, Jsc=%.4f mA/cm2, Voc=%.1f mV, FF=%.2f%%, PCE=%.2f%%\n', ...
    code,diag.nAccepted,diag.nAttempted,m.Jsc*1e3,m.Voc*1e3,m.FF*100,m.PCE);
end

summary=cell2table(rows,'VariableNames', ...
  {'Architecture','Name','Jsc_mAcm2','Voc_mV','FF_pct','PCE_pct', ...
   'AcceptedPoints','AttemptedPoints','MaxDeltaJabs_Acm2','MaxDeltaJrel', ...
   'MaxElectronResidual','MaxHoleResidual','MaxPoissonResidual'});
writetable(summary,fullfile(resultsDir,'stageA_comparison.csv'));
save(fullfile(resultsDir,'stageA_comparison.mat'),'summary','runs','Vgrid', ...
  'Pin','d1nm','d2nm','Rbar');

writeMarkdownSummary(summary,d1nm,d2nm,Rbar,Pin, ...
  fullfile(reportsDir,'paper_stageA_comparison.md'));
plotPaperStageAJV(runs,fullfile(reportsDir,'paper_stageA_jv_comparison.png'), ...
  fullfile(reportsDir,'paper_stageA_jv_comparison.pdf'));

disp(summary(:,1:8));
end

function [curve,diag]=strictSweep(dev,G,Vgrid)
N=dev.N;
states=cell(numel(Vgrid),1);
accepted=false(numel(Vgrid),1);
J=nan(numel(Vgrid),1);
dJabs=nan(numel(Vgrid),1); dJrel=nan(numel(Vgrid),1);
eRes=nan(numel(Vgrid),1); hRes=nan(numel(Vgrid),1); pRes=nan(numel(Vgrid),1);

% Start slightly above zero, then return to V=0 using reverse continuation.
idxStart=find(abs(Vgrid-0.02)<1e-12,1);
seed=struct('psi',dev.psiEq,'phin',zeros(N,1),'phip',zeros(N,1));
r=gummelSolveQFPAnalytic(dev,Vgrid(idxStart),G,struct( ...
  'psiInit',seed.psi,'phinInit',seed.phin,'phipInit',seed.phip));
states{idxStart}=r;
if r.isFinite
  r0=gummelSolveQFPAnalytic(dev,0,G,struct( ...
    'psiInit',r.psi,'phinInit',r.phin,'phipInit',r.phip));
  states{1}=r0;
end

prev=r;
for k=idxStart+1:numel(Vgrid)
  if ~prev.isFinite, break; end
  rk=gummelSolveQFPAnalytic(dev,Vgrid(k),G,struct( ...
    'psiInit',prev.psi,'phinInit',prev.phin,'phipInit',prev.phip));
  states{k}=rk;
  if rk.isFinite, prev=rk; end
end

for k=1:numel(Vgrid)
  if isempty(states{k}), continue; end
  r=states{k};
  accepted(k)=qfpResultAccepted(r);
  if r.isFinite
    J(k)=r.J;
    dJabs(k)=r.currentConservation.deltaJabs;
    dJrel(k)=r.currentConservation.rJ;
    eRes(k)=r.continuityResiduals.electron.maxRelGlobal;
    hRes(k)=r.continuityResiduals.hole.maxRelGlobal;
    pRes(k)=r.poissonResidual;
  end
end

curve=table(Vgrid,J,accepted,dJabs,dJrel,eRes,hRes,pRes, ...
  'VariableNames',{'V_V','J_Acm2','Accepted','DeltaJabs_Acm2','DeltaJrel', ...
  'ElectronResidual','HoleResidual','PoissonResidual'});
diag=struct('nAccepted',sum(accepted),'nAttempted',sum(~cellfun(@isempty,states)), ...
  'maxDeltaJabs',safeMax(dJabs(accepted)),'maxDeltaJrel',safeMax(dJrel(accepted)), ...
  'maxElectronResidual',safeMax(eRes(accepted)), ...
  'maxHoleResidual',safeMax(hRes(accepted)), ...
  'maxPoissonResidual',safeMax(pRes(accepted)));
end

function v=safeMax(x)
if isempty(x), v=NaN; else, v=max(x,[],'omitnan'); end
end

function [d1best,d2best,Rbar]=optimizeCommonARC(am15)
absOpt=loadGaAsNPOptical(3.0);
lam=(300:10:1000)';
irr=interp1(am15(:,1),am15(:,2),lam,'linear',0);
E=1239.842./lam;
alpha=interp1(absOpt.E,absOpt.alpha,E,'linear',0);
alpha(E<absOpt.Emin)=0;
alpha(E>absOpt.Emax)=absOpt.alpha(end);
absorptance=1-exp(-alpha*1500e-7);
weight=(irr./max(E,1e-6)).*absorptance;
d1range=60:5:160; d2range=40:5:140;
Rbar=Inf; d1best=NaN; d2best=NaN;
for d1=d1range
  for d2=d2range
    R=zeros(size(lam));
    for k=1:numel(lam)
      R(k)=tmmReflectance(lam(k), ...
        [1.0,1.38,2.35,gapRefractiveIndexReal(lam(k))],[d1,d2]);
    end
    rw=sum(weight.*R)/sum(weight);
    if rw<Rbar, Rbar=rw; d1best=d1; d2best=d2; end
  end
end
end

function makeJVFigure(runs,pngFile,pdfFile)
fig=figure('Visible','off','Color','w','Position',[100 100 1000 680]);
ax=axes(fig); hold(ax,'on');
colors=lines(numel(runs));
styles={'-','--','-.',':','-','--','-.'};
for k=1:numel(runs)
  c=runs(k).curve;
  ok=c.Accepted & isfinite(c.J_Acm2);
  plot(ax,c.V_V(ok),c.J_Acm2(ok)*1e3,'LineWidth',1.8, ...
    'Color',colors(k,:),'LineStyle',styles{k},'DisplayName',runs(k).code);
end
yline(ax,0,'k:','HandleVisibility','off');
xlabel(ax,'Voltage (V)'); ylabel(ax,'Current density (mA cm^{-2})');
title(ax,'Controlled AM1.5G comparison of GaAsPN architectures');
legend(ax,'Location','southwest','NumColumns',2);
grid(ax,'on'); box(ax,'on');
set(ax,'FontName','Arial','FontSize',11,'LineWidth',0.8);
exportgraphics(fig,pngFile,'Resolution',300);
exportgraphics(fig,pdfFile,'ContentType','vector');
close(fig);
end

function writeMarkdownSummary(T,d1,d2,Rbar,Pin,file)
fid=fopen(file,'w'); cleanup=onCleanup(@() fclose(fid));
fprintf(fid,'# Stage A controlled architecture comparison\n\n');
fprintf(fid,'Common AM1.5G input power: %.3f W/m2. Common ARC: MgF2 %.1f nm / ZnS %.1f nm (weighted R=%.5f).\n\n',Pin,d1,d2,Rbar);
fprintf(fid,'| Architecture | Jsc (mA/cm2) | Voc (mV) | FF (%%) | PCE (%%) | Accepted |\n');
fprintf(fid,'|---|---:|---:|---:|---:|---:|\n');
for k=1:height(T)
  fprintf(fid,'| %s | %.4f | %.1f | %.2f | %.2f | %d/%d |\n', ...
    T.Architecture(k),T.Jsc_mAcm2(k),T.Voc_mV(k),T.FF_pct(k), ...
    T.PCE_pct(k),T.AcceptedPoints(k),T.AttemptedPoints(k));
end
fprintf(fid,'\nOnly canonically accepted points enter the reported metrics and J-V curves; rejected points are omitted from the figure and disclosed through the accepted/attempted counts.\n');
end
