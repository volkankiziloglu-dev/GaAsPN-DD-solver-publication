function p = solveHoleContinuity(dev, psi, pLeftBC, pRightBC, Gprofile, nPrev, pPrev, Sback, p0Back)
% Sback, p0Back (OPTIONAL): Robin (surface-recombination-velocity) back BC,
% symmetric counterpart to solveElectronContinuity's Sback treatment -- needed
% when the Abs region is n-type (holes are the back-side minority carrier).
% If omitted, falls back to the original pure-Dirichlet BC (p(N)=pRightBC),
% preserving old behavior exactly for the standard p-type-absorber polarity.
c = constants();
N = dev.N; h = dev.h; muP = dev.muP; tau = dev.tau; ni = dev.ni;
% Holes are driven by Ev(x), not Ec(x) -- use deltaEv (= deltaEc + dEg), NOT
% deltaEc, for the SG drift argument. See buildDevice.m and
% MATHEMATICAL_MODEL.md Sec. A.2/I.1. (Fixed bug: this previously used
% dev.deltaEc, correct only when dEg=0, i.e. a homojunction.)
psiEff = psi + dev.psiOffsetP;
a = zeros(N,1); b = zeros(N,1); cc = zeros(N,1); d = zeros(N,1);
a(1) = 0; b(1) = 1; cc(1) = 0; d(1) = pLeftBC;
if nargin < 8 || isempty(Sback)
  a(N) = 0; b(N) = 1; cc(N) = 0; d(N) = pRightBC;
else
  hL1 = h(N-1);
  DpL1 = c.VT * 0.5*(muP(N)+muP(N-1));
  dPsiL1 = (psiEff(N)-psiEff(N-1))/c.VT;
  BL_p1 = sgBernoulli(dPsiL1); BL_m1 = sgBernoulli(-dPsiL1);
  % Derived by matching this file's own interior flux term (coefPim1 =
  % -(DpL/hL)*BL_p, coefPi's from-left part = +(DpL/hL)*BL_m) to a
  % Sback*(p(N)-p0Back) recombination sink at the back contact.
  a(N) = -(DpL1/hL1)*BL_p1;
  b(N) = (DpL1/hL1)*BL_m1 + Sback;
  cc(N) = 0;
  d(N) = Sback * p0Back;
end

idx = 2:N-1;
hL = h(idx-1); hR = h(idx);
DpL = c.VT*0.5*(muP(idx)+muP(idx-1));
DpR = c.VT*0.5*(muP(idx)+muP(idx+1));
dPsiL = (psiEff(idx)-psiEff(idx-1))/c.VT;
dPsiR = (psiEff(idx+1)-psiEff(idx))/c.VT;
hAvg = 0.5*(hL+hR);

BL_p = sgBernoulli(dPsiL); BL_m = sgBernoulli(-dPsiL);
BR_p = sgBernoulli(dPsiR); BR_m = sgBernoulli(-dPsiR);

coefPim1 = -(DpL./hL).*BL_p./hAvg;
coefPi = (DpL./hL).*BL_m./hAvg + (DpR./hR).*BR_p./hAvg;
coefPip1 = -(DpR./hR).*BR_m./hAvg;

n_i = nPrev(idx);
% R = R_SRH + R_rad (see solveElectronContinuity.m for the matching electron
% side and MATHEMATICAL_MODEL.md Section D for why this is not a SNS J02 term).
% FIX (external review, see solveElectronContinuity.m for the full
% derivation): constant term is ni^2/denom (+Brad*ni^2), not Rcoef.*ni --
% correct only in the special case n_i=ni. Sign of Rconst in d stays "+"
% (same role as the old Rcoef.*ni term played), only its VALUE changes.
[RcoefSRH,RconstSRH] = bulkSRHLinearization( ...
  dev,nPrev,pPrev,idx,'hole');
Rcoef = RcoefSRH + dev.Brad(idx).*n_i;
Rconst = RconstSRH + dev.Brad(idx).*ni(idx).^2;

% Interface recombination R_int (see solveElectronContinuity.m for the
% matching electron side and MATHEMATICAL_MODEL.md Section J).
% SIGN NOTE: unlike Rcoef.*ni (which is +ADDED in both n and p equations --
% but that sign is tied to the SPECIFIC SRH linearization identity, R =
% Rcoef*(carrier-ni), NOT a general rule), R_int is a standalone, fully
% explicit loss with NO paired implicit b-term. Because this equation's
% diagonal (coefPi+Rcoef) is POSITIVE (opposite sign convention from
% solveElectronContinuity's negative diagonal), a pure loss must be
% SUBTRACTED here to reduce p, not added -- confirmed empirically via
% test_interface_recomb_sanity.m + debug_Rint_values.m: with "+" here, p
% at the ETL-side junction node ran away (2e11 -> 3.9e16 cm^-3) as Sint
% increased, a positive-feedback blowup, not a physical recombination loss.
assertIfaceModelSupported(dev, mfilename);   % 'four'/'cross' only in assembleResidualQFP (dd_solver v1.1.1)
Rint_perNode = interfaceRecombination(dev, nPrev, pPrev, idx, hAvg);

a(idx) = coefPim1;
cc(idx) = coefPip1;
b(idx) = coefPi + Rcoef;
d(idx) = Gprofile(idx) + Rconst - Rint_perNode;

% NOT (N.16 arastirmasi -- bkz. solveElectronContinuity.m'deki ayni not
% icin tam gerekce): iki duzeltme denendi (kosegen on-kosullandirma,
% robustTridiagSolve.m'nin kismi-pivotlamali \ cozumu), IKISI DE
% basarisiz oldu -- sorun pivot-korumasi degil, ham yogunluk
% degiskenlerinin kendi kosul sayisindan kaynaklaniyor. thomasSolve.m'ye
% GERI DONULDU.
p = thomasSolve(a, b, cc, d);
end
