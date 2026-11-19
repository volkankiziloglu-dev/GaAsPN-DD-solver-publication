function summary = run_paper_indoor_production()
%RUN_PAPER_INDOOR_PRODUCTION Canonical 7x4x2 indoor architecture sweep.
% CIE LED-B1/B3/B5 and FL11 spectra are normalized directly to 100/1000
% lux. Uses the common Stage-A ARC and the declared full recombination
% scenario (SRH + radiative + Auger). Legacy FL19-labelled results are not
% read or overwritten.

assertCanonicalPath();
warning('off','MATLAB:nearlySingularMatrix');
warning('off','MATLAB:singularMatrix');
root=fileparts(fileparts(mfilename('fullpath')));
outDir=fullfile(root,'results','paper_indoor_production');
if ~isfolder(outDir), mkdir(outDir); end

arc=load(fullfile(root,'results','paper_stageA','stageA_comparison.mat'),'d1nm','d2nm');
arcFn=@(lam) tmmReflectance(lam, ...
  [1.0,1.38,2.35,gapRefractiveIndexReal(lam)],[arc.d1nm,arc.d2nm]);
codes=["HJ-0","HJ-U","HJ-G","F-HJ","RBSF-GaP","RBSF-GaNP","RHJ"];
modes=["2700K","4000K","6500K","4000F"];
luxLevels=[100,1000];
Vgrid=(0:0.02:1.56)';
Bmat=paperRadiativeCoefficients(300);

nRuns=numel(codes)*numel(modes)*numel(luxLevels);
rows=cell(nRuns,17); runs=struct([]); rr=0;
for im=1:numel(modes)
  for il=1:numel(luxLevels)
    spec=buildCIEIndoorSpectrum(modes(im),luxLevels(il));
    Pin=spectrumPower(spec);
    for ia=1:numel(codes)
      rr=rr+1;
      fprintf('\n[%d/%d] %s | %s | %d lux\n',rr,nRuns,codes(ia),modes(im),luxLevels(il));
      [dev,meta]=buildPaperArchitecture(codes(ia), ...
        BradAbs=Bmat.GaAsPN_cm3ps,BradGaP=Bmat.GaP_cm3ps, ...
        BradGaNP=Bmat.GaPN_Ege18_cm3ps, ...
        CnAugerAbs=3e-30,CpAugerAbs=3e-30, ...
        CnAugerGaP=3e-30,CpAugerGaP=3e-30, ...
        CnAugerGaNP=3e-30,CpAugerGaNP=3e-30);
      G=deviceGenAndTransmitML(dev,spec,arcFn);
      [curve,diag,gate]=indoorSweep(dev,G,Vgrid,Pin);
      if isfinite(gate.Jsc) && isfinite(gate.Voc) && isfinite(gate.FF) && isfinite(gate.PCE)
        m=struct('Jsc',gate.Jsc,'Voc',gate.Voc,'FF',gate.FF,'PCE',gate.PCE, ...
          'Pmax',gate.PCE/100*Pin/1e4);
      else
        m=struct('Jsc',NaN,'Voc',NaN,'FF',NaN,'PCE',NaN,'Pmax',NaN);
      end
      rows(rr,:)={codes(ia),meta.architectureName,modes(im),luxLevels(il),Pin, ...
        m.Jsc*1e6,m.Voc*1e3,m.FF*100,m.PCE,m.Pmax*1e3, ...
        diag.nAccepted,diag.nAttempted,diag.maxDeltaJabs,diag.maxDeltaJrel, ...
        diag.maxElectronResidual,diag.maxHoleResidual,gate.curveClass};
      runs(rr).code=codes(ia); runs(rr).mode=modes(im); runs(rr).lux=luxLevels(il); %#ok<AGROW>
      runs(rr).curve=curve; runs(rr).metrics=m; runs(rr).meta=meta;
      stem=sprintf('%s_%s_%dlux_jv.csv',codes(ia),modes(im),luxLevels(il));
      writetable(curve,fullfile(outDir,stem));
      fprintf('  %s; critical accepted %d/%d; Jsc=%.3f uA/cm2, Voc=%.1f mV, FF=%.2f%%, PCE=%.2f%%\n', ...
        gate.curveClass,diag.nAccepted,diag.nAttempted,m.Jsc*1e6,m.Voc*1e3,m.FF*100,m.PCE);
    end
  end
end

summary=cell2table(rows,'VariableNames',{'Architecture','Name','Spectrum', ...
 'Lux','Pin_Wm2','Jsc_uAcm2','Voc_mV','FF_pct','PCE_pct','Pmax_mWcm2', ...
 'AcceptedPoints','AttemptedPoints','MaxDeltaJabs_Acm2','MaxDeltaJrel', ...
 'MaxElectronResidual','MaxHoleResidual','CurveClass'});
writetable(summary,fullfile(outDir,'indoor_production_summary.csv'));
save(fullfile(outDir,'indoor_production_summary.mat'),'summary','runs','Vgrid', ...
 'arc','codes','modes','luxLevels','-v7.3');
writeReport(summary,fullfile(root,'reports','paper_indoor_production.md'));
disp(summary(:,{'Architecture','Spectrum','Lux','Jsc_uAcm2','Voc_mV','FF_pct','PCE_pct','AcceptedPoints'}));
end

function [curve,diag,gate]=indoorSweep(dev,G,Vgrid,Pin)
N=dev.N; idxStart=find(abs(Vgrid-0.02)<1e-12,1);
seed=struct('psi',dev.psiEq,'phin',zeros(N,1),'phip',zeros(N,1));
r002=gummelSolveQFPAnalytic(dev,0.02,G,struct( ...
 'psiInit',seed.psi,'phinInit',seed.phin,'phipInit',seed.phip));
% Resolve Jsc through the established reverse-continuation path.
prev=r002;
for v=[0.01,0.005,0]
 r0=gummelSolveQFPAnalytic(dev,v,G,struct( ...
  'psiInit',prev.psi,'phinInit',prev.phin,'phipInit',prev.phip));
 prev=r0;
end
states=repmat(r002,numel(Vgrid),1); states(1)=r0; states(idxStart)=r002;
prev=r002;
for k=idxStart+1:numel(Vgrid)
 rk=gummelSolveQFPAnalytic(dev,Vgrid(k),G,struct( ...
  'psiInit',prev.psi,'phinInit',prev.phin,'phipInit',prev.phip));
 states(k)=rk; if rk.isFinite, prev=rk; end
end
gate=indoorAcceptanceGate(dev,G,Vgrid,states,Pin);
J=gate.Jterm; accepted=gate.pointPass;
dJabs=nan(size(Vgrid)); dJrel=dJabs; eRes=dJabs; hRes=dJabs;
for k=1:numel(Vgrid)
 q=states(k);
 if q.isFinite
  dJabs(k)=q.currentConservation.deltaJabs; dJrel(k)=q.currentConservation.rJ;
  eRes(k)=gate.rnGlobal(k); hRes(k)=gate.rpGlobal(k);
 end
end
curve=table(Vgrid,J,accepted,gate.region,dJabs,dJrel,eRes,hRes,gate.rFB, ...
 'VariableNames',{'V_V','J_Acm2','Accepted','Region','DeltaJabs_Acm2', ...
 'DeltaJrel','ElectronResidual','HoleResidual','FrontBackRelativeError'});
critical=gate.region~="deepforward";
diag=struct('nAccepted',sum(accepted & critical),'nAttempted',sum(critical), ...
 'maxDeltaJabs',safeMax(dJabs(accepted & critical)), ...
 'maxDeltaJrel',safeMax(dJrel(accepted & critical)), ...
 'maxElectronResidual',safeMax(eRes(accepted & critical)), ...
 'maxHoleResidual',safeMax(hRes(accepted & critical)));
end

function v=safeMax(x)
if isempty(x), v=NaN; else, v=max(x,[],'omitnan'); end
end

function writeReport(T,file)
fid=fopen(file,'w'); c=onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'# Canonical indoor production sweep\n\n');
fprintf(fid,'CIE LED-B1/B3/B5 and FL11; common ARC; full SRH+radiative+Auger model.\n\n');
fprintf(fid,'| Spectrum | Lux | Architecture | Jsc (uA/cm2) | Voc (mV) | FF (%%) | PCE (%%) | Accepted |\n');
fprintf(fid,'|---|---:|---|---:|---:|---:|---:|---:|\n');
for k=1:height(T)
 fprintf(fid,'| %s | %d | %s | %.3f | %.1f | %.2f | %.2f | %d/%d |\n', ...
  T.Spectrum(k),T.Lux(k),T.Architecture(k),T.Jsc_uAcm2(k),T.Voc_mV(k), ...
  T.FF_pct(k),T.PCE_pct(k),T.AcceptedPoints(k),T.AttemptedPoints(k));
end
end
