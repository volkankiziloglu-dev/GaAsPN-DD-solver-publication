function J = assemblePhipJacobianAnalyticQFP(dev, psi, phin, phip, V, Gprofile)
%ASSEMBLEPHIPJACOBIANANALYTICQFP Exact d(F_phip)/d(phip), psi/phin frozen.
% Mirror of assemblePhinJacobianAnalyticQFP.m for the hole block. Sign
% conventions differ from the electron case in three places (all taken
% directly from assembleResidualQFP.m, not re-derived from scratch by
% symmetry, to avoid propagating a sign slip): (1) dp/dphip = +p/VT (vs
% dn/dphin = -n/VT); (2) JpSeg(i) = Q*DpSeg/h*(p(i)*B(+dPsiP) -
% p(i+1)*B(-dPsiP)) -- p(i) pairs with B(+.), p(i+1) with B(-.), the
% OPPOSITE assignment from JnSeg's n(i+1)/B(+.), n(i)/B(-.); (3)
% Fphip(idxC) = JpSeg(i)-JpSeg(i-1) + Q*hAvgC*(Rtot-G), a PLUS sign
% where Fphin has MINUS. Same structural ceiling note as the electron
% file applies: does not capture the phip-vs-phin cross-carrier coupling
% in interfaceRecombination.m (N.16.16 remains open).
c = constants();
N = dev.N;
h = dev.h;
[n, p] = carrierDensitiesQFP(dev, psi, phin, phip);

argP = -(psi + (dev.chi-dev.chiRef) + dev.Eg - phip) / c.VT;
EXP_CLAMP = 700;
clampedP = abs(argP) >= EXP_CLAMP;
dp_dphip = p / c.VT; dp_dphip(clampedP) = 0; % dp(i)/dphip(i), local only

idxJ = (1:N-1)';
psiEffP = psi + dev.psiOffsetP;
dPsiP = (psiEffP(idxJ+1) - psiEffP(idxJ)) / c.VT;
DpSeg = c.VT*0.5*(dev.muP(idxJ)+dev.muP(idxJ+1));
Bplus = sgBernoulli(dPsiP);   % multiplies p(i)
Bminus = sgBernoulli(-dPsiP); % multiplies p(i+1)
coefSeg = c.Q*DpSeg./h(idxJ);

% JpSeg(i) = coefSeg(i) * (p(i)*Bplus(i) - p(i+1)*Bminus(i))
dJseg_dphipI   = coefSeg .* Bplus  .* dp_dphip(idxJ);      % d/dphip(i)
dJseg_dphipIp1 = -coefSeg .* Bminus .* dp_dphip(idxJ+1);   % d/dphip(i+1)

lower = zeros(N,1); diagv = zeros(N,1); upper = zeros(N,1);
diagv(1) = 1; % Fphip(1) = phip(1)
if isfield(dev,'Sfront') && ~isempty(dev.Sfront) && isfinite(dev.Sfront) && dev.Na(1) <= 0
  % Fphip(1) = -JpSeg(1) - Q*Sfront*(p(1)-pLeft)   (optional finite front velocity, minority holes)
  diagv(1) = -dJseg_dphipI(1) - c.Q*dev.Sfront*dp_dphip(1);
  upper(1) = -dJseg_dphipIp1(1);
end

idxC = (2:N-1)';
hLc = h(idxC-1); hRc = h(idxC); hAvgC = 0.5*(hLc+hRc);

% Fphip(i) = JpSeg(i) - JpSeg(i-1) + Q*hAvgC(i)*(Rtot(i)-G(i))
lower(idxC) = -dJseg_dphipI(idxC-1);           % from -JpSeg(i-1), d/dphip(i-1)
diagv(idxC) = dJseg_dphipI(idxC) - dJseg_dphipIp1(idxC-1); % JpSeg(i) d/dphip(i) - JpSeg(i-1) d/dphip(i)
upper(idxC) = dJseg_dphipIp1(idxC);            % from JpSeg(i), d/dphip(i+1)

dRbulk_dp = bulkSRHDerivP(dev, n, p, idxC);
nn=n(idxC); pp=p(idxC); npExcess=nn.*pp-dev.ni(idxC).^2;
dRauger_dp = dev.CpAuger(idxC).*npExcess + ...
  (dev.CnAuger(idxC).*nn+dev.CpAuger(idxC).*pp).*nn;
dRtot_dp_local = dRbulk_dp + dev.Brad(idxC).*nn + dRauger_dp;
diagv(idxC) = diagv(idxC) + c.Q*hAvgC.*dRtot_dp_local.*dp_dphip(idxC); % PLUS sign (Fphip's own sign)

if dev.Na(N) > 0
  diagv(N) = 1; % Fphip(N) = phip(N) - V
else
  % Fphip(N) = JpSeg(N-1) - Q*Sback*(p(N)-pRight)
  lower(N) = dJseg_dphipI(N-1);
  diagv(N) = dJseg_dphipIp1(N-1) - c.Q*dev.Sback*dp_dphip(N);
end

rows = (1:N)';
J = sparse(rows, rows, diagv, N, N);
J = J + sparse((2:N)', (1:N-1)', lower(2:N), N, N);
J = J + sparse((1:N-1)', (2:N)', upper(1:N-1), N, N);

J = J + interfaceRecombJacobianAddP(dev, n, p, idxC, hAvgC, dp_dphip, N, c.Q);
end

function dR = bulkSRHDerivP(dev, n, p, idx)
idx = idx(:);
nn = n(idx); pp = p(idx); ni = dev.ni(idx);
if isfield(dev,'srhTauN') && ~isempty(dev.srhTauN)
  tauN = dev.srhTauN(idx,:); tauP = dev.srhTauP(idx,:);
  n1 = dev.srhN1(idx,:); p1 = dev.srhP1(idx,:);
else
  tauN = dev.tau(idx); tauP = tauN; n1 = ni; p1 = ni;
end
D = tauP.*(nn+n1) + tauN.*(pp+p1) + 1e-30;
num = nn.*pp - ni.^2;
dR = sum((nn.*D - num.*tauN)./D.^2, 2);
end

function Jadd = interfaceRecombJacobianAddP(dev, n, p, idxC, hAvgC, dp_dphip, N, Q) %#ok<INUSL>
Jadd = sparse(N,N);
ni = dev.ni;
pairs = getJunctionPairsLocal(dev);
for jk = 1:size(pairs,1)
  iE = pairs(jk,1); iA = pairs(jk,2); Sn = pairs(jk,3); Sp = pairs(jk,4);
  if Sn <= 0 || Sp <= 0, continue; end
  if isfield(dev,'ifaceModel') && any(strcmp(dev.ifaceModel,{'onesided','four','cross'}))
    if isempty(find(idxC==iE,1)) || isempty(find(idxC==iA,1)), continue; end
    if isfield(dev,'ifaceEt') && ~isempty(dev.ifaceEt), EtJ = dev.ifaceEt; else, EtJ = 0; end
    switch dev.ifaceModel
      case 'onesided'   % R = R(n_k, p_k) of node k only; whole hole sink on k
        k = dev.ifaceNodeA(jk); [n1, p1] = ifaceN1P1(dev.ifaceN0A(jk), dev.ifaceP0A(jk), EtJ); nA = n(k); pA = p(k);
        Dd = (nA + n1)/Sp + (pA + p1)/Sn; num = nA*pA - dev.ifaceN0A(jk)*dev.ifaceP0A(jk);
        Jadd(k,k) = Jadd(k,k) + Q*(nA*Dd - num/Sn)/Dd^2*dp_dphip(k);
      case 'four'       % hole sinks sE = R*pE/Np on iE, sA = R*pA/Np on iA
        [n1, p1] = ifaceN1P1(dev.ifaceN0(jk), dev.ifaceP0(jk), EtJ); pE = p(iE); pA = p(iA);
        Nn = n(iE) + n(iA); Np = pE + pA;
        Dd = (Nn + n1)/Sp + (Np + p1)/Sn; num = Nn*Np - dev.ifaceN0(jk)*dev.ifaceP0(jk); R = num/Dd;
        dR = (Nn*Dd - num/Sn)/Dd^2; wE = pE/Np; wA = pA/Np;
        Jadd(iE,iE) = Jadd(iE,iE) + Q*(dR*wE + R*pA/Np^2)*dp_dphip(iE);
        Jadd(iE,iA) = Jadd(iE,iA) + Q*(dR*wE - R*pE/Np^2)*dp_dphip(iA);
        Jadd(iA,iE) = Jadd(iA,iE) + Q*(dR*wA - R*pA/Np^2)*dp_dphip(iE);
        Jadd(iA,iA) = Jadd(iA,iA) + Q*(dR*wA + R*pE/Np^2)*dp_dphip(iA);
      case 'cross'      % R depends on p(iA); full hole sink on iA
        np0 = dev.ifaceNP0(jk); n1 = sqrt(np0); nE = n(iE); pA = p(iA);
        Dd = (nE + n1)/Sp + (pA + n1)/Sn;
        Jadd(iA,iA) = Jadd(iA,iA) + Q*(nE*Dd - (nE*pA - np0)/Sn)/Dd^2*dp_dphip(iA);
    end
    continue;
  end
  % N.19 fix: geometric mean (see assemblePhinJacobianAnalyticQFP.m's
  % matching comment and MATHEMATICAL_MODEL.md N.19).
  [wE, wA] = intRecWeights(dev, iE, iA);
  if wE == 0.5
    n_int = sqrt(n(iE)*n(iA)); p_int = sqrt(p(iE)*p(iA)); ni_int = sqrt(ni(iE)*ni(iA));
  else
    n_int = n(iE)^wE*n(iA)^wA; p_int = p(iE)^wE*p(iA)^wA; ni_int = ni(iE)^wE*ni(iA)^wA;
  end
  denom = (n_int+ni_int)/Sp + (p_int+ni_int)/Sn;
  dRint_dpint = (n_int*denom - (n_int*p_int-ni_int^2)/Sn) / denom^2;
  dpint_dpE = wE*p_int/p(iE); dpint_dpA = wA*p_int/p(iA);
  dRint_dpE = dRint_dpint * dpint_dpE * dp_dphip(iE);
  dRint_dpA = dRint_dpint * dpint_dpA * dp_dphip(iA);
  if isempty(find(idxC==iE,1)) || isempty(find(idxC==iA,1)), continue; end
  % Fphip's own sign is PLUS (opposite of Fphin's MINUS).
  Jadd(iE,iE) = Jadd(iE,iE) + Q*0.5*dRint_dpE;
  Jadd(iE,iA) = Jadd(iE,iA) + Q*0.5*dRint_dpA;
  Jadd(iA,iE) = Jadd(iA,iE) + Q*0.5*dRint_dpE;
  Jadd(iA,iA) = Jadd(iA,iA) + Q*0.5*dRint_dpA;
end
end

function pairs = getJunctionPairsLocal(dev)
pairs = zeros(0,4);
if isfield(dev,'junctionNodes') && ~isempty(dev.junctionNodes)
  for jk = 1:numel(dev.junctionNodes)
    iE = dev.junctionNodes(jk);
    pairs(end+1,:) = [iE, iE+1, dev.Sn_int_j(jk), dev.Sp_int_j(jk)]; %#ok<AGROW>
  end
elseif isfield(dev,'Sn_int') && isfield(dev,'Sp_int') && dev.Sn_int>0 && dev.Sp_int>0 && isfield(dev,'N_ETL')
  pairs = [dev.N_ETL, dev.N_ETL+1, dev.Sn_int, dev.Sp_int];
end
end
