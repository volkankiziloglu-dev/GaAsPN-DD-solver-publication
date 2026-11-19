function J = assemblePhinJacobianAnalyticQFP(dev, psi, phin, phip, V, Gprofile)
%ASSEMBLEPHINJACOBIANANALYTICQFP Exact d(F_phin)/d(phin), psi/phip frozen.
% Tridiagonal PLUS the interfaceRecombination.m cross-node term (which
% stays inside the tridiagonal band since the two junction nodes iE,iA
% are always adjacent, iA=iE+1 -- see interfaceRecombination.m). Does
% NOT capture the phin-vs-phip cross-carrier coupling interfaceRecombination.m
% also carries (n_int AND p_int both feed R_int) -- that dependence is
% structurally invisible to a phin-only block by construction (phip is
% frozen), exactly like the existing FD Jacobian for this same block.
% This is precisely why this Jacobian does not fix N.16.16's junction
% stagnation -- it is a faithful, exact derivative of what the phin-only
% block actually solves, not a wider or different problem. Any speed
% gain here is BEFORE that structural ceiling, not past it.
%
% See assembleResidualQFP.m for the residual this differentiates and
% MATHEMATICAL_MODEL.md N.16.21 for the derivation notes.
c = constants();
N = dev.N;
h = dev.h;
[n, p] = carrierDensitiesQFP(dev, psi, phin, phip);

chiRel = dev.chi - dev.chiRef;
argN = (psi + chiRel - phin) / c.VT;
EXP_CLAMP = 700;
clampedN = abs(argN) >= EXP_CLAMP;
dn_dphin = -n / c.VT; dn_dphin(clampedN) = 0; % dn(i)/dphin(i), local only

idxJ = (1:N-1)';
psiEffN = psi + dev.psiOffsetN;
dPsiN = (psiEffN(idxJ+1) - psiEffN(idxJ)) / c.VT;
DnSeg = c.VT*0.5*(dev.muN(idxJ)+dev.muN(idxJ+1));
Bplus = sgBernoulli(dPsiN);   % B(dPsiN(i))  -- multiplies n(i+1)
Bminus = sgBernoulli(-dPsiN); % B(-dPsiN(i)) -- multiplies n(i)
coefSeg = c.Q*DnSeg./h(idxJ);

% dJnSeg(i)/dphin(i) and dJnSeg(i)/dphin(i+1), i=1..N-1
dJseg_dphinI   = -coefSeg .* Bminus .* dn_dphin(idxJ);     % note dn_dphin already carries the -n/VT sign
dJseg_dphinIp1 =  coefSeg .* Bplus  .* dn_dphin(idxJ+1);

lower = zeros(N,1); diagv = zeros(N,1); upper = zeros(N,1);
diagv(1) = 1; % Fphin(1) = phin(1)
if isfield(dev,'Sfront') && ~isempty(dev.Sfront) && isfinite(dev.Sfront) && dev.Nd(1) <= 0
  % Fphin(1) = -JnSeg(1) + Q*Sfront*(n(1)-nLeft)   (optional finite front velocity, minority electrons)
  diagv(1) = -dJseg_dphinI(1) + c.Q*dev.Sfront*dn_dphin(1);
  upper(1) = -dJseg_dphinIp1(1);
end

idxC = (2:N-1)';
hLc = h(idxC-1); hRc = h(idxC); hAvgC = 0.5*(hLc+hRc);

% Flux part: Fphin(i) = JnSeg(i) - JnSeg(i-1) - Q*hAvgC(i)*(Rtot(i)-G(i))
% JnSeg(i) index into dJseg_* uses local position i (1-based into idxJ, i.e. array index i itself since idxJ=1:N-1)
lower(idxC) = -dJseg_dphinI(idxC-1);          % from -JnSeg(i-1), d/dphin(i-1)
diagv(idxC) = dJseg_dphinI(idxC) - dJseg_dphinIp1(idxC-1); % JnSeg(i) d/dphin(i) - JnSeg(i-1) d/dphin(i)
upper(idxC) = dJseg_dphinIp1(idxC);           % from JnSeg(i), d/dphin(i+1)

% Recombination part: -Q*hAvgC(i)*dRtot(i)/dphin(i)
dRbulk_dn = bulkSRHDerivN(dev, n, p, idxC);
nn=n(idxC); pp=p(idxC); npExcess=nn.*pp-dev.ni(idxC).^2;
dRauger_dn = dev.CnAuger(idxC).*npExcess + ...
  (dev.CnAuger(idxC).*nn+dev.CpAuger(idxC).*pp).*pp;
dRtot_dn_local = dRbulk_dn + dev.Brad(idxC).*pp + dRauger_dn;
diagv(idxC) = diagv(idxC) - c.Q*hAvgC.*dRtot_dn_local.*dn_dphin(idxC);

% Back boundary (i=N)
if dev.Nd(N) > 0
  diagv(N) = 1; % Fphin(N) = phin(N) - V
else
  % Fphin(N) = JnSeg(N-1) + Q*Sback*(n(N)-nRight)
  lower(N) = dJseg_dphinI(N-1);
  diagv(N) = dJseg_dphinIp1(N-1) + c.Q*dev.Sback*dn_dphin(N);
end

rows = (1:N)';
J = sparse(rows, rows, diagv, N, N);
J = J + sparse((2:N)', (1:N-1)', lower(2:N), N, N);
J = J + sparse((1:N-1)', (2:N)', upper(1:N-1), N, N);

% Interface-recombination cross-node contribution, added directly via
% sparse accumulation at the exact (iE,iA) node pair (see function below).
J = J + interfaceRecombJacobianAddN(dev, n, p, idxC, hAvgC, dn_dphin, N, c.Q);
end

function dR = bulkSRHDerivN(dev, n, p, idx)
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
dR = sum((pp.*D - num.*tauP)./D.^2, 2);
end

function [posIdx, valDiag, valOff, offTarget] = interfaceRecombDerivN(dev, n, p, idxC, hAvgC, dn_dphin) %#ok<INUSD,STOUT>
% Placeholder retained for readability of the main function's control
% flow; the actual accumulation is done by interfaceRecombJacobianAddN
% below (needs to add into a fully-assembled sparse J, not a dense
% diag/lower/upper triple, since node pairs may not be simply adjacent
% in the (idxC-1)-relative indexing used above).
posIdx = []; valDiag = []; valOff = []; offTarget = [];
end

function Jadd = interfaceRecombJacobianAddN(dev, n, p, idxC, hAvgC, dn_dphin, N, Q)
Jadd = sparse(N,N);
ni = dev.ni;
pairs = getJunctionPairs(dev);
for jk = 1:size(pairs,1)
  iE = pairs(jk,1); iA = pairs(jk,2); Sn = pairs(jk,3); Sp = pairs(jk,4);
  if Sn <= 0 || Sp <= 0, continue; end
  if isfield(dev,'ifaceModel') && any(strcmp(dev.ifaceModel,{'onesided','four','cross'}))
    if isempty(find(idxC==iE,1)) || isempty(find(idxC==iA,1)), continue; end
    if isfield(dev,'ifaceEt') && ~isempty(dev.ifaceEt), EtJ = dev.ifaceEt; else, EtJ = 0; end
    switch dev.ifaceModel
      case 'onesided'   % R = R(n_k, p_k) of node k only; whole electron sink on k
        k = dev.ifaceNodeA(jk); [n1, p1] = ifaceN1P1(dev.ifaceN0A(jk), dev.ifaceP0A(jk), EtJ); nA = n(k); pA = p(k);
        Dd = (nA + n1)/Sp + (pA + p1)/Sn; num = nA*pA - dev.ifaceN0A(jk)*dev.ifaceP0A(jk);
        Jadd(k,k) = Jadd(k,k) - Q*(pA*Dd - num/Sp)/Dd^2*dn_dphin(k);
      case 'four'       % electron sinks sE = R*nE/Nn on iE, sA = R*nA/Nn on iA
        [n1, p1] = ifaceN1P1(dev.ifaceN0(jk), dev.ifaceP0(jk), EtJ); nE = n(iE); nA = n(iA);
        Nn = nE + nA; Np = p(iE) + p(iA);
        Dd = (Nn + n1)/Sp + (Np + p1)/Sn; num = Nn*Np - dev.ifaceN0(jk)*dev.ifaceP0(jk); R = num/Dd;
        dR = (Np*Dd - num/Sp)/Dd^2; wE = nE/Nn; wA = nA/Nn;
        Jadd(iE,iE) = Jadd(iE,iE) - Q*(dR*wE + R*nA/Nn^2)*dn_dphin(iE);
        Jadd(iE,iA) = Jadd(iE,iA) - Q*(dR*wE - R*nE/Nn^2)*dn_dphin(iA);
        Jadd(iA,iE) = Jadd(iA,iE) - Q*(dR*wA - R*nA/Nn^2)*dn_dphin(iE);
        Jadd(iA,iA) = Jadd(iA,iA) - Q*(dR*wA + R*nE/Nn^2)*dn_dphin(iA);
      case 'cross'      % R depends on n(iE); full electron sink on iE
        np0 = dev.ifaceNP0(jk); n1 = sqrt(np0); nE = n(iE); pA = p(iA);
        Dd = (nE + n1)/Sp + (pA + n1)/Sn;
        Jadd(iE,iE) = Jadd(iE,iE) - Q*(pA*Dd - (nE*pA - np0)/Sp)/Dd^2*dn_dphin(iE);
    end
    continue;
  end
  % N.19 fix: n_int/p_int use the GEOMETRIC mean (matches
  % interfaceRecombination.m's own N.19 fix -- arithmetic mean was a
  % real bug, gave nonzero R_int even at true equilibrium; see
  % MATHEMATICAL_MODEL.md N.19 for the algebraic proof).
  [wE, wA] = intRecWeights(dev, iE, iA);
  if wE == 0.5
    n_int = sqrt(n(iE)*n(iA)); p_int = sqrt(p(iE)*p(iA)); ni_int = sqrt(ni(iE)*ni(iA));
  else
    n_int = n(iE)^wE*n(iA)^wA; p_int = p(iE)^wE*p(iA)^wA; ni_int = ni(iE)^wE*ni(iA)^wA;
  end
  denom = (n_int+ni_int)/Sp + (p_int+ni_int)/Sn;
  % dRint/dn_int = [p_int*denom - (n_int*p_int-ni_int^2)/Sp] / denom^2
  dRint_dnint = (p_int*denom - (n_int*p_int-ni_int^2)/Sp) / denom^2;
  % Geometric-mean chain rule: d(n_int)/d(n(iE)) = n_int/(2*n(iE)),
  % d(n_int)/d(n(iA)) = n_int/(2*n(iA)) (replaces the old constant-0.5
  % arithmetic-mean derivative).
  dnint_dnE = wE*n_int/n(iE); dnint_dnA = wA*n_int/n(iA);
  dRint_dnE = dRint_dnint * dnint_dnE * dn_dphin(iE);
  dRint_dnA = dRint_dnint * dnint_dnA * dn_dphin(iA);
  % R_int_perNode(posE) += 0.5*Rint/hAvg(posE); posA similarly.
  % Fphin(iE) -= Q*hAvgC(at iE)*Rtot(iE), and Rtot includes +0.5*Rint/hAvg(iE)
  % so the hAvg cancels exactly: contribution to Fphin(iE) is -Q*0.5*Rint,
  % independent of hAvg (matches interfaceRecombination.m's own conservative
  % per-node-hAvg-normalized-then-hAvg-weighted-back construction).
  hIdxE = find(idxC==iE,1); hIdxA = find(idxC==iA,1);
  if isempty(hIdxE) || isempty(hIdxA), continue; end
  Jadd(iE,iE) = Jadd(iE,iE) - Q*0.5*dRint_dnE;
  Jadd(iE,iA) = Jadd(iE,iA) - Q*0.5*dRint_dnA;
  Jadd(iA,iE) = Jadd(iA,iE) - Q*0.5*dRint_dnE;
  Jadd(iA,iA) = Jadd(iA,iA) - Q*0.5*dRint_dnA;
end
end

function pairs = getJunctionPairs(dev)
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
