function summary=run_paper_voltage_step_convergence()
%RUN_PAPER_VOLTAGE_STEP_CONVERGENCE Compare 0.02 V and 0.01 V grids.
assertCanonicalPath(); warning('off','MATLAB:nearlySingularMatrix'); warning('off','MATLAB:singularMatrix');
root=fileparts(fileparts(mfilename('fullpath'))); out=fullfile(root,'results','paper_voltage_step');
if ~isfolder(out), mkdir(out); end
codes=["HJ-0","HJ-U","HJ-G","F-HJ","RBSF-GaP","RBSF-GaNP","RHJ"];
steps=[0.02,0.01]; am15=loadRealAM15G(); Pin=spectrumPower(am15);
a=load(fullfile(root,'results','paper_stageA','stageA_comparison.mat'),'d1nm','d2nm');
arcFn=@(lam)tmmReflectance(lam,[1,1.38,2.35,gapRefractiveIndexReal(lam)],[a.d1nm,a.d2nm]);
rows=cell(14,9); r=0;
for i=1:numel(codes)
 for j=1:numel(steps)
  r=r+1; fprintf('[%d/14] %s dV=%.3f V\n',r,codes(i),steps(j));
  dev=fullDevice(codes(i)); G=deviceGenAndTransmitML(dev,am15,arcFn);
  V=(0:steps(j):1.70)'; [curve,diag,m]=runAcceptedJVSweep(dev,G,V,Pin);
  rows(r,:)={codes(i),steps(j),m.Jsc*1e3,m.Voc*1e3,m.FF*100,m.PCE, ...
    diag.nAccepted,diag.nAttempted,all(curve.Accepted)};
  writetable(curve,fullfile(out,sprintf('%s_dV_%g.csv',codes(i),steps(j))));
 end
end
summary=cell2table(rows,'VariableNames',{'Architecture','VoltageStep_V','Jsc_mAcm2', ...
 'Voc_mV','FF_pct','PCE_pct','AcceptedPoints','AttemptedPoints','AllAccepted'});
summary.dJsc_pct=nan(height(summary),1); summary.dVoc_mV=nan(height(summary),1);
summary.dFF_pp=nan(height(summary),1); summary.dPCE_pp=nan(height(summary),1);
summary.Pass=false(height(summary),1);
for i=1:numel(codes)
 i2=find(summary.Architecture==codes(i)&summary.VoltageStep_V==0.02);
 i1=find(summary.Architecture==codes(i)&summary.VoltageStep_V==0.01);
 d=[100*(summary.Jsc_mAcm2(i2)-summary.Jsc_mAcm2(i1))/summary.Jsc_mAcm2(i1), ...
  summary.Voc_mV(i2)-summary.Voc_mV(i1),summary.FF_pct(i2)-summary.FF_pct(i1), ...
  summary.PCE_pct(i2)-summary.PCE_pct(i1)];
 summary{i2,{'dJsc_pct','dVoc_mV','dFF_pp','dPCE_pp'}}=d;
 summary.Pass(i2)=abs(d(1))<0.5&&abs(d(2))<5&&abs(d(3))<0.5&&abs(d(4))<0.10;
end
writetable(summary,fullfile(out,'voltage_step_convergence.csv'));
save(fullfile(out,'voltage_step_convergence.mat'),'summary'); disp(summary);
end

function dev=fullDevice(code)
Bmat=paperRadiativeCoefficients(300);
dev=buildPaperArchitecture(code,BradAbs=Bmat.GaAsPN_cm3ps, ...
 BradGaP=Bmat.GaP_cm3ps,BradGaNP=Bmat.GaPN_Ege18_cm3ps, ...
 CnAugerAbs=3e-30,CpAugerAbs=3e-30,CnAugerGaP=3e-30,CpAugerGaP=3e-30, ...
 CnAugerGaNP=3e-30,CpAugerGaNP=3e-30);
end
