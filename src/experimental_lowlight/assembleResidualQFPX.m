function [F, n, p, JnSeg, JpSeg] = assembleResidualQFPX(dev, psi, phin, phip, V, Gprofile)
%ASSEMBLERESIDUALQFPX DENEYSEL (experimental_lowlight) kopya: src/assembleResidualQFP.m ile satir
% satir ayni; YALNIZ Scharfetter-Gummel segment akilari dev.fluxModeX ile secilir:
%   'prod'   (varsayilan) : uretim, n(i+1)B(x)-n(i)B(-x), sgBernoulli (bit-ozdes)
%   'expm1'               : ayni fark formu, sgBernoulliExpm1 cekirdegi
%   'excess'              : fazla-tasiyici (denge-cikarilmis) form. n(i+1)/n(i)=exp(x-dphin/VT) ve
%                           B(x)/B(-x)=exp(-x) ozdesliklerinden denge bileseni analitik olarak
%                           cikarilir: Jn = qD/h*n(i)*B(-x)*expm1(-dphin/VT),
%                                      Jp = qD/h*p(i+1)*B(-xp)*expm1(-dphip/VT).
%                           Karanlik dengede (dphi=0) akim TAM sifir; buyuk surukleme/difuzyon
%                           terimlerinin farki hic alinmaz. (Heteroarayuzde Nc/Nv orani
%                           dev.psiOffsetN/P icinde sadelesir; bkz. src/sgFluxQFPStable.m.)
% Uretim kodu degistirilmedi.
c = constants();
N = dev.N;
h = dev.h;

[n, p] = carrierDensitiesQFP(dev, psi, phin, phip);

idxP = (2:N-1)';
epsL = c.EPS0*0.5*(dev.eps(idxP)+dev.eps(idxP-1));
epsR = c.EPS0*0.5*(dev.eps(idxP)+dev.eps(idxP+1));
hL = h(idxP-1); hR = h(idxP); hAvgP = 0.5*(hL+hR);
rho = c.Q*(p(idxP)-n(idxP)+dev.Nd(idxP)-dev.Na(idxP));

Fpsi = zeros(N,1);
Fpsi(idxP) = (epsR.*(psi(idxP+1)-psi(idxP))./hR - epsL.*(psi(idxP)-psi(idxP-1))./hL)./hAvgP + rho;
Fpsi(1) = psi(1) - dev.psiEq(1);
Fpsi(N) = psi(N) - (dev.psiEq(N) + V);

idxJ = (1:N-1)';
psiEffN = psi + dev.psiOffsetN;
psiEffP = psi + dev.psiOffsetP;
dPsiN = (psiEffN(idxJ+1)-psiEffN(idxJ))/c.VT;
dPsiP = (psiEffP(idxJ+1)-psiEffP(idxJ))/c.VT;
DnSeg = c.VT*0.5*(dev.muN(idxJ)+dev.muN(idxJ+1));
DpSeg = c.VT*0.5*(dev.muP(idxJ)+dev.muP(idxJ+1));
mode = 'prod'; if isfield(dev,'fluxModeX') && ~isempty(dev.fluxModeX), mode = dev.fluxModeX; end
switch mode
  case 'prod'
    JnSeg = c.Q*DnSeg./h(idxJ) .* (n(idxJ+1).*sgBernoulli(dPsiN) - n(idxJ).*sgBernoulli(-dPsiN));
    JpSeg = c.Q*DpSeg./h(idxJ) .* (p(idxJ).*sgBernoulli(dPsiP) - p(idxJ+1).*sgBernoulli(-dPsiP));
  case 'expm1'
    JnSeg = c.Q*DnSeg./h(idxJ) .* (n(idxJ+1).*sgBernoulliExpm1(dPsiN) - n(idxJ).*sgBernoulliExpm1(-dPsiN));
    JpSeg = c.Q*DpSeg./h(idxJ) .* (p(idxJ).*sgBernoulliExpm1(dPsiP) - p(idxJ+1).*sgBernoulliExpm1(-dPsiP));
  case 'excess'
    dphin = phin(idxJ+1)-phin(idxJ); dphip = phip(idxJ+1)-phip(idxJ);
    JnSeg = c.Q*DnSeg./h(idxJ) .* n(idxJ).*sgBernoulliExpm1(-dPsiN).*expm1(-dphin/c.VT);
    JpSeg = c.Q*DpSeg./h(idxJ) .* p(idxJ+1).*sgBernoulliExpm1(-dPsiP).*expm1(-dphip/c.VT);
  case 'excessB0'
    % ayni fazla-tasiyici formu, uretim Bernoulli cekirdegi (sgBernoulli) ile
    dphin = phin(idxJ+1)-phin(idxJ); dphip = phip(idxJ+1)-phip(idxJ);
    JnSeg = c.Q*DnSeg./h(idxJ) .* n(idxJ).*sgBernoulli(-dPsiN).*expm1(-dphin/c.VT);
    JpSeg = c.Q*DpSeg./h(idxJ) .* p(idxJ+1).*sgBernoulli(-dPsiP).*expm1(-dphip/c.VT);
  otherwise
    error('assembleResidualQFPX:mode','bilinmeyen fluxModeX: %s', mode);
end

idxC = (2:N-1)';
hLc = h(idxC-1); hRc = h(idxC); hAvgC = 0.5*(hLc+hRc);
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

Fphin = zeros(N,1);
Fphip = zeros(N,1);
Fphin(idxC) = (JnSeg(idxC) - JnSeg(idxC-1)) - c.Q*hAvgC.*(RtotN - Gprofile(idxC));
Fphip(idxC) = (JpSeg(idxC) - JpSeg(idxC-1)) + c.Q*hAvgC.*(RtotP - Gprofile(idxC));

Sfront = Inf; if isfield(dev,'Sfront') && ~isempty(dev.Sfront), Sfront = dev.Sfront; end
Fphin(1) = phin(1);
Fphip(1) = phip(1);
if isfinite(Sfront)
  if dev.Na(1) > 0, pLeft = dev.Na(1); else, pLeft = dev.ni(1)^2/dev.Nd(1); end
  nLeft = dev.ni(1)^2/pLeft;
  if dev.Nd(1) <= 0, Fphin(1) = -JnSeg(1) + c.Q*Sfront*(n(1)-nLeft); end
  if dev.Na(1) <= 0, Fphip(1) = -JpSeg(1) - c.Q*Sfront*(p(1)-pLeft); end
end

Sback = dev.Sback;
if dev.Na(N) > 0, pRight = dev.Na(N); else, pRight = dev.ni(N)^2/dev.Nd(N); end
nRight = dev.ni(N)^2/pRight;
if dev.Nd(N) > 0
  Fphin(N) = phin(N) - V;
else
  Fphin(N) = JnSeg(N-1) + c.Q*Sback*(n(N)-nRight);
end
if dev.Na(N) > 0
  Fphip(N) = phip(N) - V;
else
  Fphip(N) = JpSeg(N-1) - c.Q*Sback*(p(N)-pRight);
end

F = zeros(3*N,1);
F(1:3:end) = Fpsi;
F(2:3:end) = Fphin;
F(3:3:end) = Fphip;
end
