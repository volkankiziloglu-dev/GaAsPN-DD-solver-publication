function n = solveElectronContinuity(dev, psi, nLeftBC, nRightBC, Gprofile, nPrev, pPrev, Sback, n0Back)
% Scharfetter-Gummel discretization (vectorized). Front: ohmic Dirichlet (majority).
% Back: Robin BC via Sback (surface recombination velocity) IF Sback/n0Back
% are supplied; otherwise Dirichlet (n(N)=nRightBC), mirroring
% solveHoleContinuity.m's optional-Sback pattern. FIX (external review,
% Madde 6): previously Sback/n0Back were REQUIRED, so electrons always got
% Robin treatment at the back even when they are the MAJORITY carrier there
% (n-type absorber) -- asymmetric with how holes are already dispatched
% conditionally in gummelSolve.m/gummelSolveML.m. Existing 2-region callers
% (common/gummelSolve.m) still always pass Sback explicitly for the
% standard p-type-absorber polarity (electrons minority there, Robin
% correct), so this is exact backward compatibility, not a behavior change
% for any existing validated script.
c = constants();
N = dev.N; h = dev.h; muN = dev.muN; tau = dev.tau; ni = dev.ni;
psiEff = psi + dev.psiOffsetN;
a = zeros(N,1); b = zeros(N,1); cc = zeros(N,1); d = zeros(N,1);
a(1) = 0; b(1) = 1; cc(1) = 0; d(1) = nLeftBC;

if nargin < 8 || isempty(Sback)
  a(N) = 0; b(N) = 1; cc(N) = 0; d(N) = nRightBC;
else
  hL1 = h(N-1);
  DnL1 = c.VT * 0.5*(muN(N)+muN(N-1));
  dPsiL1 = (psiEff(N)-psiEff(N-1))/c.VT;
  BL_p1 = sgBernoulli(dPsiL1); BL_m1 = sgBernoulli(-dPsiL1);
  % At the right boundary electron particle loss is opposite to the
  % conventional electron-current direction:
  %   Jn/q = -Sback*(n-n0).
  % Moving the surface term to the matrix gives +Sback on the diagonal
  % and +Sback*n0 on the right-hand side.  The former -S/-S*n0 signs
  % crossed the transport diagonal at large S and produced negative n.
  a(N) = -(DnL1/hL1)*BL_m1;
  b(N) = (DnL1/hL1)*BL_p1 + Sback;
  cc(N) = 0;
  d(N) = Sback * n0Back;
end

idx = 2:N-1;
hL = h(idx-1); hR = h(idx);
DnL = c.VT*0.5*(muN(idx)+muN(idx-1));
DnR = c.VT*0.5*(muN(idx)+muN(idx+1));
dPsiL = (psiEff(idx)-psiEff(idx-1))/c.VT;
dPsiR = (psiEff(idx+1)-psiEff(idx))/c.VT;
hAvg = 0.5*(hL+hR);

BL_p = sgBernoulli(dPsiL); BL_m = sgBernoulli(-dPsiL);
BR_p = sgBernoulli(dPsiR); BR_m = sgBernoulli(-dPsiR);

coefNim1 = (DnL./hL).*BL_m./hAvg;
coefNi = -(DnL./hL).*BL_p./hAvg - (DnR./hR).*BR_m./hAvg;
coefNip1 = (DnR./hR).*BR_p./hAvg;

p_i = pPrev(idx);
% R = R_SRH + R_rad, both linearized (p frozen at previous iterate) to
% R(n) ~= Rcoef*n - Rconst. dev.Brad=0 (default) preserves old SRH-only
% behavior exactly. See buildDevice.m and MATHEMATICAL_MODEL.md Section D
% -- NOT a Sah-Noyce-Shockley J02 term (that would double-count
% recombination already integrated volumetrically here).
% FIX (external review): the constant term is ni^2/denom (+ Brad*ni^2), NOT
% Rcoef*ni. R_SRH(n) = (n*p_prev-ni^2)/denom = (p_prev/denom)*n - ni^2/denom
% -- the two forms coincide only in the special case p_prev=ni; in doped
% regions (p_prev>>ni or <<ni) the old "Rcoef*ni" constant was wrong by a
% factor of order p_prev/ni, so the Gummel loop was converging to a
% different, physically incorrect equation, not just a less-accurate
% linearization of the right one.
[RcoefSRH,RconstSRH] = bulkSRHLinearization( ...
  dev,nPrev,pPrev,idx,'electron');
Rcoef = RcoefSRH + dev.Brad(idx).*p_i;
Rconst = RconstSRH + dev.Brad(idx).*ni(idx).^2;

% Interface recombination R_int at the ETL/Abs junction (2D/areal, separate
% from bulk SRH/radiative above) -- see interfaceRecombination.m and
% MATHEMATICAL_MODEL.md Section J. Zero vector when Sn_int/Sp_int unset
% (default), exactly recovering the "ideal interface" behavior used
% throughout this model until now.
% SIGN: "+" here because this equation's diagonal (coefNi-Rcoef) is
% NEGATIVE -- for n~d/b with b<0, a pure loss must make d LESS negative
% (i.e. ADD a positive quantity) to reduce n. solveHoleContinuity.m's
% diagonal is POSITIVE, so the same loss must be SUBTRACTED there instead
% -- see that file's note (confirmed empirically after a sign error caused
% p to run away at the junction node as Sint increased).
assertIfaceModelSupported(dev, mfilename);   % 'four'/'cross' only in assembleResidualQFP (dd_solver v1.1.1)
Rint_perNode = interfaceRecombination(dev, nPrev, pPrev, idx, hAvg);

a(idx) = coefNim1;
cc(idx) = coefNip1;
b(idx) = coefNi - Rcoef;
% Steady-state electron continuity:
%   div(Jn)/q = R - G
% with R = Rcoef*n - Rconst + Rint.  Since Rcoef*n is moved to the
% left-hand side, BOTH generation and the recombination constant enter
% the right-hand side with a minus sign.
d(idx) = -Gprofile(idx) - Rconst + Rint_perNode;

% NOT (N.16 arastirmasi): item1 dark-current dogrulamasi sirasinda, bu
% sistemin yuksek-V/guclu-uretim kosullarinda n(x)'i hacmin buyuk
% kisminda yanlislikla 1e-30 tabanina cakip akim surekliligini
% (Jtot(x)=sabit) bozdugu bulundu (bkz. MATHEMATICAL_MODEL.md N.16). IKI
% duzeltme denendi ve HER IKISI DE basarisiz oldu, GERI ALINDI:
% (1) kosegen on-kosullandirma -- tek bir dondurulmus psi anlik
%     goruntusunde calisti, TAM Gummel dongusune gomulunce calismadi;
% (2) robustTridiagSolve.m (kismi-pivotlamali yogun \ cozumu) -- hedef
%     durumda HICBIR iyilesme saglamadi VE onceden calisan baska bir
%     cihazda (2-katmanli literatur karsilastirmasi) YENI bir yakinsama
%     basarisizligi yaratti. Bu, sorunun pivot-korumasi degil, HAM
%     yogunluk degiskenlerinin (~30 mertebe dinamik aralik) COZUM
%     UZAYININ KENDI kosul sayisindan kaynaklandigini kanitliyor -- hicbir
%     dogrusal-cebir algoritmasi (pivotlu/pivotsuz) cift-duyarlikli
%     aritmetigin kendi hassasiyet sinirini asamaz. Gercek duzeltme,
%     degiskenlerin kendisini yeniden formule etmeyi (log-yogunluk/
%     Slotboom, tam-baglasik Newton icinde) gerektiriyor -- bu, ayri, cok
%     daha buyuk bir gorev. thomasSolve.m'ye GERI DONULDU.
n = thomasSolve(a, b, cc, d);
end
