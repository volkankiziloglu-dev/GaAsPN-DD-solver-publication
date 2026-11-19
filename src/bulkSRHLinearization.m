function [Rcoef,Rconst] = bulkSRHLinearization(dev,nPrev,pPrev,idx,carrier)
%BULKSRHLINEARIZATION Frozen-denominator multi-center SRH coefficients.
%
% For center j:
% Rj = (n*p-ni^2) /
%      (tauP_j*(n+n1_j) + tauN_j*(p+p1_j)).
%
% The denominator and the opposite carrier are frozen at the current
% Gummel iterate. The returned coefficients satisfy
%   electron solve: R ~= Rcoef*n - Rconst
%   hole solve:     R ~= Rcoef*p - Rconst.
%
% If the device has no explicit multi-center arrays, one mid-gap center
% with tauN=tauP=dev.tau is used. This exactly recovers the legacy model.

idx=idx(:);
n=nPrev(idx); p=pPrev(idx); ni=dev.ni(idx);
if isfield(dev,'srhTauN') && ~isempty(dev.srhTauN)
  tauN=dev.srhTauN(idx,:);
  tauP=dev.srhTauP(idx,:);
  n1=dev.srhN1(idx,:);
  p1=dev.srhP1(idx,:);
else
  tauN=dev.tau(idx); tauP=tauN;
  n1=ni; p1=ni;
end

denom=tauP.*(n+n1)+tauN.*(p+p1)+1e-30;
switch lower(string(carrier))
  case "electron"
    Rcoef=sum(p./denom,2);
  case "hole"
    Rcoef=sum(n./denom,2);
  otherwise
    error('bulkSRHLinearization:Carrier', ...
      'carrier must be electron or hole.');
end
Rconst=sum(ni.^2./denom,2);
end
