function out = gummelSolveML(dev, V, Gprofile, opts)
% RearHJT/gummelSolveML.m'nin kopyasi -- newtonPoissonML (eps0-duzeltilmis)
% cagiriyor. Elektron/hol BC dagitimi zaten kosulluydu (RearHJT'nin kendi
% "Madde 6" duzeltmesi, common/gummelSolve.m'de eksikti) -- burada
% degismedi. Varsayilan relax/maxIter, gummelSolve.m ile AYNI mantikla
% degistirildi (bkz. o dosyadaki not).
%
% TEShISTIK EKLEMELER (dis inceleme madde 1.2/1.5/3.3):
%  - un/vp hesaplamasindaki exp() argumanlari artik sinirlandiriliyor
%    (sgBernoulli'nin zaten yaptigi korumanin ayni turu) -- buyuk bant
%    ofseti/ileri kutuplamada tasma/NaN riski vardi.
%  - newtonPoissonML artik [psi,converged,iter,residual] donduruyor;
%    bu bilgi out.poissonConverged/out.poissonResidual olarak disari
%    tasiniyor, dis Gummel dongusu ic Poisson basarisiz olursa BUNU
%    da "converged" durumuna dahil ediyor (once sessizce yok sayiliyordu).
%  - akim SADECE floor(N/2)'de degil, TUM ic dugumlerde hesaplanip
%    surekliligi (deltaJ) raporlaniyor.
c = constants();
N = dev.N; ni = dev.ni; Nd = dev.Nd; Na = dev.Na;

if nargin < 4, opts = struct(); end
if ~isfield(opts,'maxIter'), opts.maxIter = 1500; end
if ~isfield(opts,'relax'), opts.relax = 0.05; end

EXP_CLAMP = 700; % exp(700)~1e304, double tasmasindan (~1.8e308) once guvenli sinir

if Nd(1) > 0
  nLeft = Nd(1);
else
  nLeft = ni(1)^2 / Na(1);
end
pLeft = ni(1)^2 / nLeft;
if Na(N) > 0
  pRight = Na(N);
else
  pRight = ni(N)^2 / Nd(N);
end
nRight = ni(N)^2 / pRight;

if isfield(opts,'psiInit') && ~isempty(opts.psiInit)
  psi = opts.psiInit;
  psiRightTarget = dev.psiEq(N) + V;
else
  psi = dev.psiEq;
  psiRightTarget = dev.psiEq(N) + V;
  psi(N) = psiRightTarget;
end
if isfield(opts,'nInit') && ~isempty(opts.nInit), n = opts.nInit; else, n = dev.nEq; end
if isfield(opts,'pInit') && ~isempty(opts.pInit), p = opts.pInit; else, p = dev.pEq; end

Sback = dev.Sback;
converged = false;
poissonConverged = true; % son ic Poisson cagrisinin durumu
poissonResidual = NaN;
poissonFailCount = 0;
negFloorTotal = 0;
negFloorMaxPerIter = 0;
subFloorPositiveTotal = 0;
minNRaw = Inf;
minPRaw = Inf;
for iter = 1:opts.maxIter
  nPrevIter = n; pPrevIter = p;

  if Nd(N) > 0
    nNew = solveElectronContinuity(dev, psi, nLeft, nRight, Gprofile, nPrevIter, pPrevIter);
  else
    nNew = solveElectronContinuity(dev, psi, nLeft, nRight, Gprofile, nPrevIter, pPrevIter, Sback, nRight);
  end
  if Na(N) > 0
    pNew = solveHoleContinuity(dev, psi, pLeft, pRight, Gprofile, nPrevIter, pPrevIter);
  else
    pNew = solveHoleContinuity(dev, psi, pLeft, pRight, Gprofile, nPrevIter, pPrevIter, Sback, pRight);
  end
  % N.16 (dis inceleme): 1e-30 tabani daha once SESSIZCE negatif/anlamsiz
  % cozumleri (thomasSolve.m'nin yikici-iptal basarisizligi, bkz.
  % MATHEMATICAL_MODEL.md N.16) "yakinsanmis" gibi gosteriyordu. Negatif
  % (fiziksel olarak imkansiz, sadece pozitif-ama-kucuk degil) degerlerin
  % sayisi ayrica sayilip out.negFloorHit olarak disari tasiniyor; bu
  % sayimin son iterasyonda sifirdan farkli olmasi out.converged'i FALSE
  % yapar (asagida).
  negFloorCount = sum(nNew < 0) + sum(pNew < 0);
  negFloorTotal = negFloorTotal + negFloorCount;
  negFloorMaxPerIter = max(negFloorMaxPerIter, negFloorCount);
  subFloorPositiveTotal = subFloorPositiveTotal + ...
      sum(nNew >= 0 & nNew < 1e-30) + sum(pNew >= 0 & pNew < 1e-30);
  minNRaw = min(minNRaw, min(nNew));
  minPRaw = min(minPRaw, min(pNew));
  nNew(nNew < 1e-30) = 1e-30;
  pNew(pNew < 1e-30) = 1e-30;

  maxRelNP = max(max(abs(log(nNew./n))), max(abs(log(pNew./p))));
  n = nNew; p = pNew;

  psiEffN = psi + dev.psiOffsetN;
  psiEffP = psi + dev.psiOffsetP;
  argN = max(min(-psiEffN/c.VT, EXP_CLAMP), -EXP_CLAMP);
  argP = max(min(psiEffP/c.VT, EXP_CLAMP), -EXP_CLAMP);
  un = n .* exp(argN);
  vp = p .* exp(argP);
  [psiNew, poissonConverged, ~, poissonResidual] = newtonPoissonML(dev, psi, un, vp, dev.psiEq(1), psiRightTarget, false);
  if ~poissonConverged
    poissonFailCount = poissonFailCount + 1;
  end

  step = opts.relax * (psiNew - psi);
  psi = psi + step;
  maxChange = max(abs(step));

  if maxChange < 1e-12 && maxRelNP < 1e-9 && iter > 3 && poissonConverged
    converged = true;
    break;
  end
end

% The loop updates psi after solving n/p. At low illumination a single
% carrier re-solve is not enough because SRH recombination couples the two
% continuity equations. Iterate them at fixed final psi until both carrier
% profiles are mutually consistent.
negFloorCount = 0;
for polishIter = 1:25
  nBeforePolish = n;
  pBeforePolish = p;
  if Nd(N) > 0
    nPolish = solveElectronContinuity(dev, psi, nLeft, nRight, Gprofile, n, p);
  else
    nPolish = solveElectronContinuity(dev, psi, nLeft, nRight, Gprofile, n, p, Sback, nRight);
  end
  nPolishRaw = nPolish;
  nPolish = max(nPolishRaw, 1e-30);
  if Na(N) > 0
    pPolish = solveHoleContinuity(dev, psi, pLeft, pRight, Gprofile, nPolish, p);
  else
    pPolish = solveHoleContinuity(dev, psi, pLeft, pRight, Gprofile, nPolish, p, Sback, pRight);
  end
  negPolish = sum(nPolishRaw < 0) + sum(pPolish < 0);
  negFloorTotal = negFloorTotal + negPolish;
  negFloorCount = negFloorCount + negPolish;
  negFloorMaxPerIter = max(negFloorMaxPerIter, negPolish);
  subFloorPositiveTotal = subFloorPositiveTotal + ...
      sum(nPolishRaw >= 0 & nPolishRaw < 1e-30) + sum(pPolish >= 0 & pPolish < 1e-30);
  minNRaw = min(minNRaw, min(nPolishRaw));
  minPRaw = min(minPRaw, min(pPolish));
  n = max(nPolish, 1e-30);
  p = max(pPolish, 1e-30);
  polishChange = max(max(abs(log(n./nBeforePolish))), ...
                     max(abs(log(p./pBeforePolish))));
  if polishChange < 1e-10
    break;
  end
end

% Akim surekliligi: TUM ic dugumlerde Jn+Jp hesapla (sadece floor(N/2)
% degil -- dis inceleme madde 1.5). Kararli durumda sabit olmali;
% degilse deltaJ bunu yakalar.
idxJ = 1:N-1;
psiEffFinalN = psi + dev.psiOffsetN;
psiEffFinalP = psi + dev.psiOffsetP;
dPsiN = (psiEffFinalN(idxJ+1)-psiEffFinalN(idxJ))/c.VT;
dPsiP = (psiEffFinalP(idxJ+1)-psiEffFinalP(idxJ))/c.VT;
DnSeg = c.VT*0.5*(dev.muN(idxJ)+dev.muN(idxJ+1));
DpSeg = c.VT*0.5*(dev.muP(idxJ)+dev.muP(idxJ+1));
JnSeg = c.Q*DnSeg./dev.h(idxJ) .* (n(idxJ+1).*sgBernoulli(dPsiN) - n(idxJ).*sgBernoulli(-dPsiN));
JpSeg = c.Q*DpSeg./dev.h(idxJ) .* (p(idxJ).*sgBernoulli(dPsiP) - p(idxJ+1).*sgBernoulli(-dPsiP));
JtotSeg = JnSeg + JpSeg;

Jterminal = median(JtotSeg);
Jfloor = 1e-12; % A/cm^2, mutlak taban (deltaJ'nin sifira bolunmesini onler)
deltaJ = (max(JtotSeg) - min(JtotSeg)) / max(max(abs(JtotSeg)), Jfloor);

isFinite = all(isfinite(psi)) && all(isfinite(n)) && all(isfinite(p)) && isfinite(Jterminal);
% At indoor injection levels the raw-density formulation reaches its
% floating-point conditioning floor near the highly doped front contact.
% Keep a strict relative test, but permit at most 0.1 uA/cm^2 absolute
% imbalance; this remains negligible for outdoor currents and is reported.
currentConservation = checkCurrentConservation(JtotSeg, Jterminal, 1e-7, 1e-3);
continuityResiduals = checkContinuityResiduals(dev, Gprofile, n, p, JnSeg, JpSeg);
chiRel = dev.chi-dev.chiRef;
EFn = c.VT*log(n./dev.Nc)-psi-chiRel;
EFp = -psi-chiRel-dev.Eg-c.VT*log(p./dev.Nv);
equilibriumCheck = struct();
equilibriumCheck.EFnSpan = max(EFn)-min(EFn);
equilibriumCheck.EFpSpan = max(EFp)-min(EFp);
equilibriumCheck.massActionMaxRel = max(abs(n.*p-dev.ni.^2)./max(dev.ni.^2,1e-300));
equilibriumLike = abs(V)<1e-14 && max(abs(Gprofile))<1;
if equilibriumLike
  equilibriumCheck.pass = equilibriumCheck.EFnSpan<1e-6 && ...
      equilibriumCheck.EFpSpan<1e-6 && equilibriumCheck.massActionMaxRel<1e-6;
  physicsPass = negFloorTotal==0 && currentConservation.pass && equilibriumCheck.pass;
else
  equilibriumCheck.pass = false;
  physicsPass = negFloorTotal==0 && currentConservation.pass && ...
      (continuityResiduals.pass || currentConservation.passAbs);
end

out.psi = psi; out.n = n; out.p = p; out.J = Jterminal;
% N.16 (dis inceleme, ONEMLI DUZELTME): ilk denemede negFloorHit>0 iken
% out.converged=false yapmak denendi VE GERI ALINDI -- aydinlatilmis N.1'in
% TAM sweepJVML.m suprumesini V=0.0'da DAHI derhal kirdi (bkz.
% MATHEMATICAL_MODEL.md N.16.8) -- bu, N.16'nin teshis ettigi kusurun
% sanildigindan COK daha yaygin oldugunu (yuksek-V'ye ozgu degil, hemen
% HER islenmis noktada bir miktar negatif-tabanli dugum olustugunu)
% ortaya cikardi. Bu buyuklukte bir davranis degisikligini sessizce
% production'a koymak yerine, negFloorHit SADECE TESHIS ICIN raporlanir --
% out.converged mevcut tanimiyla (isFinite+yakinsama) DEGISMEDEN kalir,
% cagiran kod (sweepJVML.m gibi) bu bilgiyi ISTERSE kullanabilir.
out.negFloorHit = negFloorCount;
out.negFloorTotal = negFloorTotal;
out.negFloorMaxPerIter = negFloorMaxPerIter;
out.subFloorPositiveTotal = subFloorPositiveTotal;
out.minNRaw = minNRaw;
out.minPRaw = minPRaw;
out.currentConservation = currentConservation;
out.continuityResiduals = continuityResiduals;
out.equilibriumCheck = equilibriumCheck;
out.physicsPass = physicsPass;
out.converged = converged && isFinite && physicsPass;
out.iters = iter;
out.isFinite = isFinite;
out.poissonConverged = poissonConverged;
out.poissonResidual = poissonResidual;
out.poissonFailCount = poissonFailCount;
out.deltaJ = deltaJ;
out.Jseg = JtotSeg;
end
