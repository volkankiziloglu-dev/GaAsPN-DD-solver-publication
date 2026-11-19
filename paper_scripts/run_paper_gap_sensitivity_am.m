function summary=run_paper_gap_sensitivity_am()
%RUN_PAPER_GAP_SENSITIVITY_AM GaP chi x lifetime grid with the unified extractor (all finite points); output results/paper_gap_sensitivity_am.
assertCanonicalPath(); warning('off','MATLAB:nearlySingularMatrix'); warning('off','MATLAB:singularMatrix');
root=fileparts(fileparts(mfilename('fullpath'))); out=fullfile(root,'results','paper_gap_sensitivity_am');
if ~isfolder(out), mkdir(out); end
codes=["HJ-0","HJ-U","HJ-G","F-HJ","RBSF-GaP","RBSF-GaNP","RHJ"];
chis=[3.7,3.8,4.0,4.1,4.3]; taus=[1e-12,1e-11,1e-10,1e-9,1e-8];
am15=loadRealAM15G(); Pin=spectrumPower(am15); V=(0:0.02:1.70)';
a=load(fullfile(root,'results','paper_stageA','stageA_comparison.mat'),'d1nm','d2nm');
arcFn=@(lam)tmmReflectance(lam,[1,1.38,2.35,gapRefractiveIndexReal(lam)],[a.d1nm,a.d2nm]);
Bmat=paperRadiativeCoefficients(300);
n=numel(codes)*numel(chis)*numel(taus); rows=cell(n,12); r=0;
for ic=1:numel(chis)
 for it=1:numel(taus)
  for ia=1:numel(codes)
   r=r+1; fprintf('[%d/%d] %s chi=%.2f eV tau=%.0e s\n',r,n,codes(ia),chis(ic),taus(it));
   dev=buildPaperArchitecture(codes(ia),GaPChi_eV=chis(ic),GaPTau_s=taus(it), ...
    BradAbs=Bmat.GaAsPN_cm3ps,BradGaP=Bmat.GaP_cm3ps, ...
    BradGaNP=Bmat.GaPN_Ege18_cm3ps, ...
    CnAugerAbs=3e-30,CpAugerAbs=3e-30,CnAugerGaP=3e-30,CpAugerGaP=3e-30, ...
    CnAugerGaNP=3e-30,CpAugerGaNP=3e-30);
   G=deviceGenAndTransmitML(dev,am15,arcFn); [~,diag,m]=runAcceptedJVSweep(dev,G,V,Pin,true);
   rows(r,:)={codes(ia),chis(ic),taus(it),m.Jsc*1e3,m.Voc*1e3,m.FF*100,m.PCE, ...
    diag.nAccepted,diag.nAttempted,diag.maxDeltaJrel,diag.maxElectronResidual,diag.maxHoleResidual};
   partial=cell2table(rows(1:r,:),'VariableNames',names());
   writetable(partial,fullfile(out,'gap_sensitivity_checkpoint.csv'));
  end
 end
end
summary=cell2table(rows,'VariableNames',names());
writetable(summary,fullfile(out,'gap_sensitivity_summary.csv'));
save(fullfile(out,'gap_sensitivity_summary.mat'),'summary','codes','chis','taus');
plotMaps(summary,codes,chis,taus,out); disp(summary);
end

function n=names()
n={'Architecture','GaPChi_eV','GaPTau_s','Jsc_mAcm2','Voc_mV','FF_pct','PCE_pct', ...
 'AcceptedPoints','AttemptedPoints','MaxDeltaJrel','MaxElectronResidual','MaxHoleResidual'};
end

function plotMaps(T,codes,chis,taus,out)
vars={'Jsc_mAcm2','Voc_mV','FF_pct','PCE_pct'};
labels={'J_{sc} (mA cm^{-2})','V_{oc} (mV)','FF (%)','PCE (%)'};
for iv=1:numel(vars)
 f=figure('Visible','off','Color','w','Position',[50 50 1450 720]);
 tl=tiledlayout(f,2,4,'TileSpacing','compact','Padding','compact');
 for ia=1:numel(codes)
  ax=nexttile(tl); Z=nan(numel(taus),numel(chis));
  for i=1:numel(chis), for j=1:numel(taus)
   q=T(T.Architecture==codes(ia)&T.GaPChi_eV==chis(i)&T.GaPTau_s==taus(j),:);
   if ~isempty(q), Z(j,i)=q.(vars{iv})(1); end
  end, end
  imagesc(ax,1:numel(chis),1:numel(taus),Z); axis(ax,'xy'); colorbar(ax);
  xticks(ax,1:numel(chis)); xticklabels(ax,string(chis));
  yticks(ax,1:numel(taus)); yticklabels(ax,{'10^{-12}','10^{-11}','10^{-10}','10^{-9}','10^{-8}'});
  title(ax,codes(ia)); xlabel(ax,'\chi_{GaP} (eV)'); ylabel(ax,'\tau_{GaP} (s)');
 end
 title(tl,['GaP parameter sensitivity: ' labels{iv}]);
 exportgraphics(f,fullfile(out,['gap_sensitivity_' vars{iv} '.png']),'Resolution',300);
 exportgraphics(f,fullfile(out,['gap_sensitivity_' vars{iv} '.pdf']),'ContentType','vector'); close(f);
end
end
