function summary = run_paper_recombination_scenarios()
%RUN_PAPER_RECOMBINATION_SCENARIOS Four common-physics Stage-A scenarios.
% The same proxy B/C coefficients are applied to every layer so that this
% run tests ranking robustness without introducing architecture-specific
% fitting. Material-specific coefficients remain a later uncertainty task.

assertCanonicalPath();
warning('off','MATLAB:nearlySingularMatrix');
warning('off','MATLAB:singularMatrix');
root=fileparts(fileparts(mfilename('fullpath')));
outDir=fullfile(root,'results','paper_recombination_scenarios');
reportDir=fullfile(root,'reports');
if ~exist(outDir,'dir'), mkdir(outDir); end

codes=["HJ-0","HJ-U","HJ-G","F-HJ","RBSF-GaP","RBSF-GaNP","RHJ"];
scenarioNames=["SRH","SRH+Rad","SRH+Auger","SRH+Rad+Auger"];
Bmat=paperRadiativeCoefficients(300);
Bnom=Bmat.GaAsPN_cm3ps; Cnom=3e-30;
Bvals=[0,Bnom,0,Bnom]; Cvals=[0,0,Cnom,Cnom];

am15=loadRealAM15G(); Pin=spectrumPower(am15);
[d1nm,d2nm,Rbar]=commonARC(am15);
arcFn=@(lam) tmmReflectance(lam, ...
  [1.0,1.38,2.35,gapRefractiveIndexReal(lam)],[d1nm,d2nm]);
Vgrid=(0:0.02:1.70)';
rows=cell(numel(codes)*numel(scenarioNames),11); row=0;
runs=struct('scenario',{},'architecture',{},'curve',{},'metrics',{},'diagnostics',{});

for s=1:numel(scenarioNames)
  fprintf('\n######## %s: B=%.3e, Cn=Cp=%.3e ########\n', ...
    scenarioNames(s),Bvals(s),Cvals(s));
  for a=1:numel(codes)
    row=row+1; code=codes(a);
    useRad=Bvals(s)>0;
    [dev,meta]=buildPaperArchitecture(code, ...
      BradAbs=useRad*Bmat.GaAsPN_cm3ps, ...
      BradGaP=useRad*Bmat.GaP_cm3ps, ...
      BradGaNP=useRad*Bmat.GaPN_Ege18_cm3ps, ...
      CnAugerAbs=Cvals(s),CpAugerAbs=Cvals(s), ...
      CnAugerGaP=Cvals(s),CpAugerGaP=Cvals(s), ...
      CnAugerGaNP=Cvals(s),CpAugerGaNP=Cvals(s));
    G=deviceGenAndTransmitML(dev,am15,arcFn);
    [curve,diag]=strictSweepLocal(dev,G,Vgrid);
    accepted=curve(curve.Accepted,["V_V","J_Acm2"]);
    if height(accepted)>=4 && any(abs(accepted.V_V)<1e-12) && any(accepted.J_Acm2<0)
      m=pvMetrics(table2array(accepted),Pin);
    else
      m=struct('Jsc',NaN,'Voc',NaN,'FF',NaN,'PCE',NaN,'Pmax',NaN);
    end
    runs(row)=struct('scenario',scenarioNames(s),'architecture',code, ...
      'curve',curve,'metrics',m,'diagnostics',diag); %#ok<AGROW>
    rows(row,:)={scenarioNames(s),code,Bvals(s),Cvals(s),m.Jsc*1e3, ...
      m.Voc*1e3,m.FF*100,m.PCE,diag.nAccepted,diag.nAttempted,meta.architectureName};
    writetable(curve,fullfile(outDir,scenarioNames(s)+"_"+code+"_jv.csv"));
    fprintf('%-10s accepted %d/%d, Voc=%.1f mV, PCE=%.3f%%\n', ...
      code,diag.nAccepted,diag.nAttempted,m.Voc*1e3,m.PCE);
  end
end

summary=cell2table(rows,'VariableNames',{'Scenario','Architecture','Brad_cm3ps', ...
  'CnEqCp_cm6ps','Jsc_mAcm2','Voc_mV','FF_pct','PCE_pct','AcceptedPoints', ...
  'AttemptedPoints','Name'});

base=summary(summary.Scenario=="SRH",{'Architecture','PCE_pct','Voc_mV'});
base.Properties.VariableNames={'Architecture','PCE_SRH_pct','Voc_SRH_mV'};
summary=outerjoin(summary,base,'Keys','Architecture','MergeKeys',true,'Type','left');
summary.DeltaPCE_vs_SRH_pp=summary.PCE_pct-summary.PCE_SRH_pct;
summary.DeltaVoc_vs_SRH_mV=summary.Voc_mV-summary.Voc_SRH_mV;
writetable(summary,fullfile(outDir,'recombination_scenario_comparison.csv'));
save(fullfile(outDir,'recombination_scenario_comparison.mat'),'summary','runs', ...
  'Vgrid','Pin','d1nm','d2nm','Rbar','Bnom','Bmat','Cnom');
writeReport(summary,Bnom,Cnom,Pin,d1nm,d2nm,Rbar, ...
  fullfile(reportDir,'paper_recombination_scenarios.md'));
plotScenarioPCE(summary,fullfile(reportDir,'paper_recombination_scenarios.png'), ...
  fullfile(reportDir,'paper_recombination_scenarios.pdf'));
disp(summary(:,{'Scenario','Architecture','Jsc_mAcm2','Voc_mV','FF_pct', ...
  'PCE_pct','DeltaPCE_vs_SRH_pp','AcceptedPoints','AttemptedPoints'}));
end

function [curve,diag]=strictSweepLocal(dev,G,Vgrid)
N=dev.N; states=cell(numel(Vgrid),1); accepted=false(numel(Vgrid),1);
J=nan(numel(Vgrid),1); dJabs=J; dJrel=J; eRes=J; hRes=J; pRes=J;
idxStart=find(abs(Vgrid-0.02)<1e-12,1);
r=gummelSolveQFPAnalytic(dev,Vgrid(idxStart),G,struct('psiInit',dev.psiEq, ...
  'phinInit',zeros(N,1),'phipInit',zeros(N,1)));
states{idxStart}=r;
if r.isFinite
  states{1}=gummelSolveQFPAnalytic(dev,0,G,struct('psiInit',r.psi, ...
    'phinInit',r.phin,'phipInit',r.phip));
end
prev=r;
for k=idxStart+1:numel(Vgrid)
  if ~prev.isFinite, break; end
  rk=gummelSolveQFPAnalytic(dev,Vgrid(k),G,struct('psiInit',prev.psi, ...
    'phinInit',prev.phin,'phipInit',prev.phip));
  states{k}=rk; if rk.isFinite, prev=rk; end
end
for k=1:numel(Vgrid)
  if isempty(states{k}), continue; end
  r=states{k}; accepted(k)=qfpResultAccepted(r);
  if r.isFinite
    J(k)=r.J; dJabs(k)=r.currentConservation.deltaJabs;
    dJrel(k)=r.currentConservation.rJ;
    eRes(k)=r.continuityResiduals.electron.maxRelGlobal;
    hRes(k)=r.continuityResiduals.hole.maxRelGlobal; pRes(k)=r.poissonResidual;
  end
end
curve=table(Vgrid,J,accepted,dJabs,dJrel,eRes,hRes,pRes,'VariableNames', ...
  {'V_V','J_Acm2','Accepted','DeltaJabs_Acm2','DeltaJrel', ...
  'ElectronResidual','HoleResidual','PoissonResidual'});
diag=struct('nAccepted',sum(accepted),'nAttempted',sum(~cellfun(@isempty,states)));
end

function [d1best,d2best,Rbar]=commonARC(am15)
absOpt=loadGaAsNPOptical(3.0); lam=(300:10:1000)';
irr=interp1(am15(:,1),am15(:,2),lam,'linear',0); E=1239.842./lam;
alpha=interp1(absOpt.E,absOpt.alpha,E,'linear',0);
alpha(E<absOpt.Emin)=0; alpha(E>absOpt.Emax)=absOpt.alpha(end);
weight=(irr./max(E,1e-6)).*(1-exp(-alpha*1500e-7));
Rbar=Inf; d1best=NaN; d2best=NaN;
for d1=60:5:160
  for d2=40:5:140
    R=zeros(size(lam));
    for k=1:numel(lam)
      R(k)=tmmReflectance(lam(k),[1,1.38,2.35,gapRefractiveIndexReal(lam(k))],[d1,d2]);
    end
    rw=sum(weight.*R)/sum(weight);
    if rw<Rbar, Rbar=rw; d1best=d1; d2best=d2; end
  end
end
end

function writeReport(T,B,C,Pin,d1,d2,Rbar,file)
fid=fopen(file,'w'); cleanup=onCleanup(@() fclose(fid));
fprintf(fid,'# Recombination-physics scenario comparison\n\n');
fprintf(fid,'Common proxy coefficients: B=%.3e cm3/s and Cn=Cp=%.3e cm6/s. ',B,C);
fprintf(fid,'These are sensitivity anchors, not composition-specific measured constants.\n\n');
fprintf(fid,'AM1.5G Pin=%.3f W/m2; ARC MgF2 %.1f nm / ZnS %.1f nm; weighted R=%.5f.\n\n',Pin,d1,d2,Rbar);
fprintf(fid,'| Scenario | Architecture | Jsc | Voc (mV) | FF (%%) | PCE (%%) | Delta PCE (pp) | Accepted |\n');
fprintf(fid,'|---|---|---:|---:|---:|---:|---:|---:|\n');
for k=1:height(T)
  fprintf(fid,'| %s | %s | %.4f | %.1f | %.2f | %.3f | %+.3f | %d/%d |\n', ...
    T.Scenario(k),T.Architecture(k),T.Jsc_mAcm2(k),T.Voc_mV(k),T.FF_pct(k), ...
    T.PCE_pct(k),T.DeltaPCE_vs_SRH_pp(k),T.AcceptedPoints(k),T.AttemptedPoints(k));
end
end

function plotScenarioPCE(T,pngFile,pdfFile)
scenarios=unique(T.Scenario,'stable'); codes=unique(T.Architecture,'stable');
Y=nan(numel(codes),numel(scenarios));
for i=1:numel(codes), for j=1:numel(scenarios)
  q=T.Architecture==codes(i) & T.Scenario==scenarios(j); Y(i,j)=T.PCE_pct(q);
end, end
fig=figure('Visible','off','Color','w','Position',[100 100 1100 650]);
b=bar(categorical(codes,codes),Y,'grouped');
ylabel('Power conversion efficiency (%)'); xlabel('Architecture');
title('Recombination-physics sensitivity under common AM1.5G conditions');
legend(b,cellstr(scenarios),'Location','southoutside','NumColumns',4);
grid on; box on; set(gca,'FontName','Arial','FontSize',11);
exportgraphics(fig,pngFile,'Resolution',300); exportgraphics(fig,pdfFile,'ContentType','vector'); close(fig);
end
