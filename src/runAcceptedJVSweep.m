function [curve,diag,m] = runAcceptedJVSweep(dev,G,Vgrid,Pin,useAllFinite)
%RUNACCEPTEDJVSWEEP Canonical AM1.5G sweep and accepted-point metrics. With useAllFinite = true the metrics use all finite points of the sweep
% (unified extractor, as for the indoor curves) instead of the points that pass the acceptance gate; the default (false) is unchanged.
if nargin < 5, useAllFinite = false; end
N=dev.N; states=cell(numel(Vgrid),1); accepted=false(numel(Vgrid),1);
J=nan(numel(Vgrid),1); dJabs=J; dJrel=J; eRes=J; hRes=J; pRes=J;
idxStart=find(abs(Vgrid-0.02)<1e-12,1);
if isempty(idxStart), error('Vgrid must contain 0.02 V for continuation.'); end
seed=struct('psi',dev.psiEq,'phin',zeros(N,1),'phip',zeros(N,1));
r=gummelSolveQFPAnalytic(dev,Vgrid(idxStart),G,struct( ...
 'psiInit',seed.psi,'phinInit',seed.phin,'phipInit',seed.phip));
states{idxStart}=r;
if r.isFinite
  prev0=r;
  rev=sort(Vgrid(Vgrid<0.02),'descend');
  for v=rev(:)'
    r0=gummelSolveQFPAnalytic(dev,v,G,struct( ...
      'psiInit',prev0.psi,'phinInit',prev0.phin,'phipInit',prev0.phip));
    states{find(abs(Vgrid-v)<1e-12,1)}=r0; prev0=r0;
  end
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
  hRes(k)=q.continuityResiduals.hole.maxRelGlobal; pRes(k)=q.poissonResidual;
 end
end
curve=table(Vgrid,J,accepted,dJabs,dJrel,eRes,hRes,pRes,'VariableNames', ...
 {'V_V','J_Acm2','Accepted','DeltaJabs_Acm2','DeltaJrel', ...
 'ElectronResidual','HoleResidual','PoissonResidual'});
if useAllFinite, ok=curve(isfinite(curve.J_Acm2),["V_V","J_Acm2"]); else, ok=curve(curve.Accepted,["V_V","J_Acm2"]); end
if height(ok)>=4 && any(abs(ok.V_V)<1e-12) && any(ok.J_Acm2<0)
 m=pvMetrics(table2array(ok),Pin);
else
 m=struct('Jsc',NaN,'Voc',NaN,'FF',NaN,'PCE',NaN,'Pmax',NaN);
end
diag=struct('nAccepted',sum(accepted),'nAttempted',sum(~cellfun(@isempty,states)), ...
 'maxDeltaJabs',safeMax(dJabs(accepted)),'maxDeltaJrel',safeMax(dJrel(accepted)), ...
 'maxElectronResidual',safeMax(eRes(accepted)),'maxHoleResidual',safeMax(hRes(accepted)), ...
 'maxPoissonResidual',safeMax(pRes(accepted)));
end

function v=safeMax(x)
if isempty(x), v=NaN; else, v=max(x,[],'omitnan'); end
end
