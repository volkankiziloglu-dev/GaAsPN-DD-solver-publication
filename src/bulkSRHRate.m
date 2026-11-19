function R = bulkSRHRate(dev,n,p,idx)
%BULKSRHRATE Evaluate single- or multi-center bulk SRH recombination.
idx=idx(:);
nn=n(idx); pp=p(idx); ni=dev.ni(idx);
if isfield(dev,'srhTauN') && ~isempty(dev.srhTauN)
  tauN=dev.srhTauN(idx,:);
  tauP=dev.srhTauP(idx,:);
  n1=dev.srhN1(idx,:);
  p1=dev.srhP1(idx,:);
else
  tauN=dev.tau(idx); tauP=tauN;
  n1=ni; p1=ni;
end
denom=tauP.*(nn+n1)+tauN.*(pp+p1)+1e-30;
R=sum((nn.*pp-ni.^2)./denom,2);
end
