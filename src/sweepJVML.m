function [pts, status] = sweepJVML(dev, Gprofile, Vmin, Vmax, Vstep, initialState)
% RearHJT/sweepJVML.m'nin kopyasi -- gummelSolveML (eps0-duzeltilmis)
% cagiriyor, ARTI V=0->Vmin isinma supurmesi (eps0 duzeltmesinden sonra
% Vmin sifirdan uzaksa tek buyuk bir Newton sicramasi thomasSolve.m'de
% katastrofik iptale/NaN'a yol acabiliyordu -- bkz. MATHEMATICAL_MODEL.md
% Bolum K madde 12/13).
%
% Monoton-J durdurma esigi ORANSAL (N.7'de bulunan mutlak-tolerans
% hatasinin duzeltmesi): JtolAbs = max(1e-7, 1e-3*J0).
%
% DUZELTME (dis inceleme madde 1.1/1.4): r.converged && r.isFinite HER
% noktada kontrol ediliyor; basarisiz olursa onceki yakinsanmis
% noktadan baslayarak gerilim adimi ikiye bolunerek (continuation/
% bisection, en fazla maxSubdiv kademe) tekrar denenir.
%
% DUZELTME (dis inceleme madde 2.4, REVIZE EDILMIS): ilk onerilen
% "deltaJ'yi kabul esigine kat" fikri DENENDI ve YANLIS ciktigi
% gosterildi -- debug_deltaJ_calibration.m/debug_deltaJ_profile_good.m
% ile, ZATEN dogrulanmis N.1 homojonksiyonunun BILE aydinlatma altinda
% deltaJ~1.4-1.8 (%140-180!) verdigi bulundu -- Config A'nin "bozuk"
% ~0.1-0.4'undan BUYUK. Neden: aydinlatma altinda Jtot(x)'in DUZ olmasi
% zaten fiziksel olarak beklenmez (Beer-Lambert soguru profili G(x)'i
% derinlikle degistirir, bu yuzden toplam akim her noktada yerel
% uretim-rekombinasyon butcesine gore degisir -- SADECE karanlikta
% duz olmasi beklenir). deltaJ bu haliyle iyi/kotu ayirt edemiyor.
% Bunun yerine, Config A'yi GERCEKTEN yakalayan sey kullaniliyor:
% fotoakim, sogurulan MAKSIMUM foton akisini (JscMax) FIZIKSEL olarak
% ASAMAZ -- saglam bir korunum siniri, keyfi bir esik degil. Bu sinir
% SADECE J>0 (fotoakim-baskin, Jsc civari) rejiminde uygulanir; V,
% Voc'u gectiginde J guclu negatif (karanlik-akim-baskin) olur ve
% |J| JscMax'i asabilir -- bu NORMAL'dir, sinirlanmaz.
%
% Tam, uretim/rekombinasyon-butcesi-farkinda bir surekluk-residual
% kontrolu (dis incelemenin asil istedigi) AYRI, daha buyuk bir gorev
% olarak ertelendi.
%
% DUZELTME (dis inceleme madde 2.6): ikinci (opsiyonel) cikis olarak
% `status` struct'i eklendi -- cagiran kod artik supurmenin NEDEN
% durdugunu (normal Vmax'a ulasti mi, yakinsamadi mi, akim fiziksel
% siniri asti mi, monoton-J ihlali mi) acikca gorebilir. Mevcut
% `pts = sweepJVML(...)` seklindeki eski cagrilar DEGISMEDEN calismaya
% devam eder (MATLAB'da fazladan cikis argumani geriye-uyumludur).
c = constants();
JscMax = c.Q * trapz(dev.x, Gprofile); % A/cm^2, J>0 rejiminde fiziksel ust sinir
JscMaxMargin = 1.05; % %5 pay -- kaba sweep/interpolasyon gurultusu icin

status = struct('completed',false,'stopReason','','lastVoltage',NaN,'nPoints',0);

warmupStep = 0.002;
if nargin < 6 || isempty(initialState)
  prev = struct();
else
  required={'psi','n','p','J'};
  if ~all(isfield(initialState,required))
    error('sweepJVML:InitialStateFields', ...
      'initialState must contain psi, n, p, and J.');
  end
  prev=initialState;
end
prevV = 0;
if abs(Vmin) > warmupStep && nargin < 6
  Vwarm = 0:(warmupStep*sign(Vmin)):Vmin;
  if Vwarm(end) ~= Vmin, Vwarm(end+1) = Vmin; end
  for k = 1:numel(Vwarm)
    V = Vwarm(k);
    if isfield(prev,'psi')
      opts = struct('psiInit',prev.psi,'nInit',prev.n,'pInit',prev.p);
    else
      opts = struct();
    end
    r = gummelSolveML(dev, V, Gprofile, opts);
    if ~r.isFinite
      error('sweepJVML: isinma supurmesinde V=%.4f''te NaN/Inf', V);
    end
    if ~r.converged
      error('sweepJVML: isinma supurmesinde V=%.4f''te yakinsamadi (iters=%d, poissonConverged=%d)', V, r.iters, r.poissonConverged);
    end
    if abs(V) < 1e-12 && abs(r.J) > JscMax*JscMaxMargin
      error('sweepJVML: isinma supurmesinde V=%.4f''te J=%.4e, fiziksel JscMax=%.4e''i asiyor', V, r.J, JscMax);
    end
    prev = r;
    prevV = V;
  end
end
if ~isfield(prev,'psi')
  r = gummelSolveML(dev, Vmin, Gprofile, struct());
  if ~r.isFinite || ~r.converged
    error('sweepJVML: V=%.4f''te ilk cozum yakinsamadi/NaN (iters=%d)', Vmin, r.iters);
  end
  prev = r; prevV = Vmin;
end

Vs = Vmin:Vstep:Vmax;
n = numel(Vs);
pts = zeros(n,2);
nValid = 0;
Jprev = Inf;
J0 = abs(prev.J);
JtolAbs = max(1e-7, 1e-3*J0);
maxSubdiv = 4;

for k = 1:n
  Vtarget = Vs(k);
  if abs(Vtarget - prevV) < 1e-12
    r = prev; ok = true; failReason = '';
  else
    [r, ok, failReason] = stepToConverged(dev, Gprofile, prev, prevV, Vtarget, maxSubdiv, JscMax, JscMaxMargin);
  end
  if ~ok
    status.stopReason = failReason;
    status.lastVoltage = prevV;
    break;
  end
  if k > 1 && r.J > Jprev + JtolAbs
    status.stopReason = 'nonmonotonic';
    status.lastVoltage = prevV;
    break;
  end
  nValid = k;
  pts(k,:) = [Vtarget, r.J];
  Jprev = r.J;
  prev = r; prevV = Vtarget;
end
pts = pts(1:nValid,:);
status.nPoints = nValid;
if isempty(status.stopReason)
  status.completed = true;
  status.stopReason = 'reached_Vmax';
  status.lastVoltage = prevV;
end
if nValid == 0
  status.stopReason = 'empty_curve';
end
end

function [r, ok, failReason] = stepToConverged(dev, Gprofile, prevState, prevV, targetV, depthLeft, JscMax, JscMaxMargin)
opts = struct('psiInit',prevState.psi,'nInit',prevState.n,'pInit',prevState.p);
r = gummelSolveML(dev, targetV, Gprofile, opts);
physicallyValid = ~(abs(targetV) < 1e-12 && abs(r.J) > JscMax*JscMaxMargin);
if r.isFinite && r.converged && physicallyValid
  ok = true; failReason = '';
  return;
end
if ~r.isFinite
  failReason = 'nonfinite';
elseif ~physicallyValid
  failReason = 'current_nonconservation';
else
  failReason = 'nonconvergence';
end
if depthLeft <= 0
  ok = false;
  return;
end
Vmid = 0.5*(prevV + targetV);
[rMid, okMid, failReason] = stepToConverged(dev, Gprofile, prevState, prevV, Vmid, depthLeft-1, JscMax, JscMaxMargin);
if ~okMid
  ok = false;
  return;
end
[r, ok, failReason] = stepToConverged(dev, Gprofile, rMid, Vmid, targetV, depthLeft-1, JscMax, JscMaxMargin);
end
