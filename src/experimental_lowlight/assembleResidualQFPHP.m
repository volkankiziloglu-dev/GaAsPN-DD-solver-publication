function [F, n, p, JnSeg, JpSeg] = assembleResidualQFPHP(dev, S, V, Gprofile, kind, nDigits)
%ASSEMBLERESIDUALQFPHP DENEYSEL genisletilmis hassasiyetli QFP artigi.
% Durum, iki-parcali (double-double benzeri) vektorlerle verilir:
%   psi = S.psi + S.psiLo, phin = S.phin + S.phinLo, phip = S.phip + S.phipLo
% kind:
%  'vpadiff'  : URETIM ayriklastirmasi (fark formu n(i+1)B(x)-n(i)B(-x)); yogunluklar, Bernoulli
%               ve akilar vpa (nDigits basamak) ile hesaplanir, sonra double'a yuvarlanir. Bu,
%               ayriklastirmayi degistirmeden yalniz aritmetik hassasiyetini artirir.
%  'ddexcess' : fazla-tasiyici formu (assembleResidualQFPX 'excess') ama dphi ve dpsi farklari
%               hi/lo parcalarindan kompanse edilerek alinir (durum nicemlemesi ~ulp(phi)/VT
%               yerine ~ulp(dphi)/VT); tamamen double, hizli.
% Rekombinasyon/uretim/Poisson terimleri uretimdeki gibi double'da (yuvarlanmis n,p ile).
if nargin < 6 || isempty(nDigits), nDigits = 40; end
c = constants();
N = dev.N; h = dev.h; idxJ = (1:N-1)';
chiRel = dev.chi - dev.chiRef;
DnSeg = c.VT*0.5*(dev.muN(idxJ)+dev.muN(idxJ+1));
DpSeg = c.VT*0.5*(dev.muP(idxJ)+dev.muP(idxJ+1));
coefN = c.Q*DnSeg./h(idxJ); coefP = c.Q*DpSeg./h(idxJ);

switch kind
  case 'vpadiff'
    d0 = digits; digits(nDigits); cleanupD = onCleanup(@() digits(d0));
    P = @(a) vpa(sym(a,'f'));
    psiS = P(S.psi) + P(S.psiLo); phinS = P(S.phin) + P(S.phinLo); phipS = P(S.phip) + P(S.phipLo);
    VT = P(c.VT); chiS = P(chiRel);
    nS = P(dev.Nc) .* exp((psiS + chiS - phinS)/VT);
    pS = P(dev.Nv) .* exp(-(psiS + chiS + P(dev.Eg) - phipS)/VT);
    psiEffN = psiS + P(dev.psiOffsetN); psiEffP = psiS + P(dev.psiOffsetP);
    xN = (psiEffN(idxJ+1)-psiEffN(idxJ))/VT; xP = (psiEffP(idxJ+1)-psiEffP(idxJ))/VT;
    BxN = bernVpa(xN); BmxN = bernVpa(-xN); BxP = bernVpa(xP); BmxP = bernVpa(-xP);
    JnSeg = double(P(coefN) .* (nS(idxJ+1).*BxN - nS(idxJ).*BmxN));
    JpSeg = double(P(coefP) .* (pS(idxJ).*BxP - pS(idxJ+1).*BmxP));
    n = double(nS); p = double(pS);
  case 'ddexcess'
    argN = ((S.psi - S.phin) + chiRel + (S.psiLo - S.phinLo))/c.VT;
    argP = -((S.psi - S.phip) + chiRel + dev.Eg + (S.psiLo - S.phipLo))/c.VT;
    n = dev.Nc.*exp(argN); p = dev.Nv.*exp(argP);
    dPsiN = ((S.psi(idxJ+1)-S.psi(idxJ)) + (dev.psiOffsetN(idxJ+1)-dev.psiOffsetN(idxJ)) + (S.psiLo(idxJ+1)-S.psiLo(idxJ)))/c.VT;
    dPsiP = ((S.psi(idxJ+1)-S.psi(idxJ)) + (dev.psiOffsetP(idxJ+1)-dev.psiOffsetP(idxJ)) + (S.psiLo(idxJ+1)-S.psiLo(idxJ)))/c.VT;
    dphin = (S.phin(idxJ+1)-S.phin(idxJ)) + (S.phinLo(idxJ+1)-S.phinLo(idxJ));
    dphip = (S.phip(idxJ+1)-S.phip(idxJ)) + (S.phipLo(idxJ+1)-S.phipLo(idxJ));
    JnSeg = coefN .* n(idxJ).*sgBernoulliExpm1(-dPsiN).*expm1(-dphin/c.VT);
    JpSeg = coefP .* p(idxJ+1).*sgBernoulliExpm1(-dPsiP).*expm1(-dphip/c.VT);
  otherwise
    error('assembleResidualQFPHP:kind','bilinmeyen kind %s', kind);
end

% --- geri kalan: uretim artiginin aynisi (double) ---
idxP = (2:N-1)';
epsL = c.EPS0*0.5*(dev.eps(idxP)+dev.eps(idxP-1));
epsR = c.EPS0*0.5*(dev.eps(idxP)+dev.eps(idxP+1));
hL = h(idxP-1); hR = h(idxP); hAvgP = 0.5*(hL+hR);
rho = c.Q*(p(idxP)-n(idxP)+dev.Nd(idxP)-dev.Na(idxP));
dR = (S.psi(idxP+1)-S.psi(idxP)) + (S.psiLo(idxP+1)-S.psiLo(idxP));
dL = (S.psi(idxP)-S.psi(idxP-1)) + (S.psiLo(idxP)-S.psiLo(idxP-1));
Fpsi = zeros(N,1);
Fpsi(idxP) = (epsR.*dR./hR - epsL.*dL./hL)./hAvgP + rho;
Fpsi(1) = (S.psi(1) - dev.psiEq(1)) + S.psiLo(1);
Fpsi(N) = (S.psi(N) - (dev.psiEq(N) + V)) + S.psiLo(N);

idxC = (2:N-1)';
hAvgC = 0.5*(h(idxC-1)+h(idxC));
npExcess = n(idxC).*p(idxC)-dev.ni(idxC).^2;
Rbulk = bulkSRHRate(dev, n, p, idxC) + dev.Brad(idxC).*npExcess + ...
  (dev.CnAuger(idxC).*n(idxC)+dev.CpAuger(idxC).*p(idxC)).*npExcess;
if isfield(dev,'ifaceModel') && any(strcmp(dev.ifaceModel,{'cross','four'}))
  RtotN = Rbulk + interfaceRecombination(dev, n, p, idxC, hAvgC, 'n');
  RtotP = Rbulk + interfaceRecombination(dev, n, p, idxC, hAvgC, 'p');
else
  Rint = interfaceRecombination(dev, n, p, idxC, hAvgC);
  RtotN = Rbulk + Rint; RtotP = RtotN;
end
Fphin = zeros(N,1); Fphip = zeros(N,1);
Fphin(idxC) = (JnSeg(idxC) - JnSeg(idxC-1)) - c.Q*hAvgC.*(RtotN - Gprofile(idxC));
Fphip(idxC) = (JpSeg(idxC) - JpSeg(idxC-1)) + c.Q*hAvgC.*(RtotP - Gprofile(idxC));
Sfront = Inf; if isfield(dev,'Sfront') && ~isempty(dev.Sfront), Sfront = dev.Sfront; end
Fphin(1) = S.phin(1) + S.phinLo(1);
Fphip(1) = S.phip(1) + S.phipLo(1);
if isfinite(Sfront)
  if dev.Na(1) > 0, pLeft = dev.Na(1); else, pLeft = dev.ni(1)^2/dev.Nd(1); end
  nLeft = dev.ni(1)^2/pLeft;
  if dev.Nd(1) <= 0, Fphin(1) = -JnSeg(1) + c.Q*Sfront*(n(1)-nLeft); end
  if dev.Na(1) <= 0, Fphip(1) = -JpSeg(1) - c.Q*Sfront*(p(1)-pLeft); end
end
Sback = dev.Sback;
if dev.Na(N) > 0, pRight = dev.Na(N); else, pRight = dev.ni(N)^2/dev.Nd(N); end
nRight = dev.ni(N)^2/pRight;
if dev.Nd(N) > 0, Fphin(N) = (S.phin(N) - V) + S.phinLo(N); else, Fphin(N) = JnSeg(N-1) + c.Q*Sback*(n(N)-nRight); end
if dev.Na(N) > 0, Fphip(N) = (S.phip(N) - V) + S.phipLo(N); else, Fphip(N) = JpSeg(N-1) - c.Q*Sback*(p(N)-pRight); end
F = zeros(3*N,1); F(1:3:end) = Fpsi; F(2:3:end) = Fphin; F(3:3:end) = Fphip;
end

function B = bernVpa(x)
% x/(exp(x)-1) in vpa; tam sifir noktalarda B=1.
z = (double(x) == 0);
xs = x; xs(z) = 1;          % gecici, 0/0 onlemek icin
B = xs./(exp(xs)-1);
B(z) = 1;
end
