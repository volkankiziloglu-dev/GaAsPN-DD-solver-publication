function summary = run_paper_mesh_convergence()
%RUN_PAPER_MESH_CONVERGENCE Discretization check for representative devices.
% HJ-U covers the homojunction family and RHJ covers the strongest band
% discontinuity. The 2x mesh is the fine reference. Acceptance limits
% compare the production 1x mesh against 2x: |dPCE|<0.10 percentage point,
% |dVoc|<5 mV, relative |dJsc|<0.5%, and |dFF|<0.5 percentage point.

assertCanonicalPath();
warning('off','MATLAB:nearlySingularMatrix');
warning('off','MATLAB:singularMatrix');
root = fileparts(fileparts(mfilename('fullpath')));
outDir = fullfile(root,'results','paper_mesh_convergence');
if ~isfolder(outDir), mkdir(outDir); end

base = load(fullfile(root,'results','paper_stageA','stageA_comparison.mat'), ...
    'd1nm','d2nm');
am15 = loadRealAM15G();
Pin = spectrumPower(am15);
arcFn = @(lam) tmmReflectance(lam, ...
    [1.0,1.38,2.35,gapRefractiveIndexReal(lam)],[base.d1nm,base.d2nm]);
Bmat = paperRadiativeCoefficients(300);

codes = ["HJ-U","RHJ"];
scales = [0.5,1,2];
Vgrid = (0:0.02:1.70)';
rows = cell(numel(codes)*numel(scales),13);
r = 0;
for i = 1:numel(codes)
  for j = 1:numel(scales)
    r = r+1;
    fprintf('\nMesh convergence: %s, scale %.1fx\n',codes(i),scales(j));
    [dev,meta] = buildPaperArchitecture(codes(i),meshScale=scales(j), ...
      BradAbs=Bmat.GaAsPN_cm3ps,BradGaP=Bmat.GaP_cm3ps, ...
      BradGaNP=Bmat.GaPN_Ege18_cm3ps, ...
      CnAugerAbs=3e-30,CpAugerAbs=3e-30, ...
      CnAugerGaP=3e-30,CpAugerGaP=3e-30, ...
      CnAugerGaNP=3e-30,CpAugerGaNP=3e-30);
    G = deviceGenAndTransmitML(dev,am15,arcFn);
    [curve,diag] = strictSweep(dev,G,Vgrid);
    accepted = curve(curve.Accepted,["V_V","J_Acm2"]);
    if height(accepted)>=4 && any(abs(accepted.V_V)<1e-12) && any(accepted.J_Acm2<0)
      m = pvMetrics(table2array(accepted),Pin);
    else
      m = struct('Jsc',NaN,'Voc',NaN,'FF',NaN,'PCE',NaN);
    end
    rows(r,:) = {codes(i),scales(j),meta.totalMeshNodes,m.Jsc*1e3, ...
      m.Voc*1e3,m.FF*100,m.PCE,diag.nAccepted,diag.nAttempted, ...
      diag.maxDeltaJabs,diag.maxDeltaJrel,diag.maxElectronResidual, ...
      diag.maxHoleResidual};
    writetable(curve,fullfile(outDir,sprintf('%s_mesh_%g.csv',codes(i),scales(j))));
  end
end

summary = cell2table(rows,'VariableNames',{'Architecture','MeshScale','Nodes', ...
  'Jsc_mAcm2','Voc_mV','FF_pct','PCE_pct','AcceptedPoints','AttemptedPoints', ...
  'MaxDeltaJabs_Acm2','MaxDeltaJrel','MaxElectronResidual','MaxHoleResidual'});

summary.dJsc_1x_to_fine_pct = nan(height(summary),1);
summary.dVoc_1x_to_fine_mV = nan(height(summary),1);
summary.dFF_1x_to_fine_pp = nan(height(summary),1);
summary.dPCE_1x_to_fine_pp = nan(height(summary),1);
summary.Pass1xVsFine = false(height(summary),1);
for i = 1:numel(codes)
  i1 = find(summary.Architecture==codes(i) & summary.MeshScale==1);
  i2 = find(summary.Architecture==codes(i) & summary.MeshScale==2);
  dJ = 100*(summary.Jsc_mAcm2(i1)-summary.Jsc_mAcm2(i2))/summary.Jsc_mAcm2(i2);
  dV = summary.Voc_mV(i1)-summary.Voc_mV(i2);
  dF = summary.FF_pct(i1)-summary.FF_pct(i2);
  dP = summary.PCE_pct(i1)-summary.PCE_pct(i2);
  summary{i1,{'dJsc_1x_to_fine_pct','dVoc_1x_to_fine_mV', ...
    'dFF_1x_to_fine_pp','dPCE_1x_to_fine_pp'}} = [dJ,dV,dF,dP];
  summary.Pass1xVsFine(i1) = abs(dJ)<0.5 && abs(dV)<5 && abs(dF)<0.5 && abs(dP)<0.10;
end
writetable(summary,fullfile(outDir,'mesh_convergence_summary.csv'));
save(fullfile(outDir,'mesh_convergence_summary.mat'),'summary','Vgrid','Pin');
disp(summary);
end

function [curve,diag]=strictSweep(dev,G,Vgrid)
N=dev.N; states=cell(numel(Vgrid),1); accepted=false(numel(Vgrid),1);
J=nan(numel(Vgrid),1); dJabs=J; dJrel=J; eRes=J; hRes=J;
idxStart=find(abs(Vgrid-0.02)<1e-12,1);
seed=struct('psi',dev.psiEq,'phin',zeros(N,1),'phip',zeros(N,1));
r=gummelSolveQFPAnalytic(dev,Vgrid(idxStart),G,struct( ...
  'psiInit',seed.psi,'phinInit',seed.phin,'phipInit',seed.phip));
states{idxStart}=r;
if r.isFinite
  states{1}=gummelSolveQFPAnalytic(dev,0,G,struct( ...
    'psiInit',r.psi,'phinInit',r.phin,'phipInit',r.phip));
end
prev=r;
for k=idxStart+1:numel(Vgrid)
  if ~prev.isFinite, break; end
  rk=gummelSolveQFPAnalytic(dev,Vgrid(k),G,struct( ...
    'psiInit',prev.psi,'phinInit',prev.phin,'phipInit',prev.phip));
  states{k}=rk; if rk.isFinite, prev=rk; end
end
for k=1:numel(Vgrid)
  if isempty(states{k}), continue; end
  q=states{k}; accepted(k)=qfpResultAccepted(q);
  if q.isFinite
    J(k)=q.J; dJabs(k)=q.currentConservation.deltaJabs;
    dJrel(k)=q.currentConservation.rJ;
    eRes(k)=q.continuityResiduals.electron.maxRelGlobal;
    hRes(k)=q.continuityResiduals.hole.maxRelGlobal;
  end
end
curve=table(Vgrid,J,accepted,dJabs,dJrel,eRes,hRes,'VariableNames', ...
 {'V_V','J_Acm2','Accepted','DeltaJabs_Acm2','DeltaJrel', ...
  'ElectronResidual','HoleResidual'});
diag=struct('nAccepted',sum(accepted),'nAttempted',sum(~cellfun(@isempty,states)), ...
 'maxDeltaJabs',safeMax(dJabs(accepted)),'maxDeltaJrel',safeMax(dJrel(accepted)), ...
 'maxElectronResidual',safeMax(eRes(accepted)),'maxHoleResidual',safeMax(hRes(accepted)));
end

function v=safeMax(x)
if isempty(x), v=NaN; else, v=max(x,[],'omitnan'); end
end
