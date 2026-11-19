function x = thomasSolve(a, b, c, d)
% Solve tridiagonal system a(i)x(i-1)+b(i)x(i)+c(i)x(i+1)=d(i). a(1),c(n) unused.
%
% DUZELTME (dis inceleme madde 1.3): ONCEDEN sifira yakin/kucuk pivot
% (m = b(i)-a(i)*cp(i-1)) icin HICBIR koruma yoktu -- Inf/NaN/buyuk iptal
% hatasi uretebilirdi (warm-up rampasinin gerekcesinde gecmiste tam da
% bu turden bir kotu-kosullanma yasandigi belirtilmisti, bkz.
% MATHEMATICAL_MODEL.md Bolum K). Pivot mutlak degeri, HER SATIRIN
% KENDI b(i) olcegine gore (satir-yerel) normlanmis bir esigin altina
% duserse UYARI verilip NaN vektoru donduruluyor.
%
% ONEMLI: ilk versiyon TEK, GLOBAL bir esik kullanmisti (norm(b,inf)
% temelli) -- bu, sinir-kosulu satirlarinda (b(1)=1, kasitli olarak
% normlanmis Dirichlet denklemi) YANLIS ALARM veriyordu, cunku ic
% dugumlerin difuzyon katsayilari (b~1e13-1e15 mertebesinde olabilir,
% bkz. K.12'nin Case B izlemesi) global normu domine edip b(1)=1'i
% "sifira yakin" gibi gosteriyordu (capraz-kontrol scripti calistirilinca
% HEMEN yakalandi -- debug_thomassolve_crosscheck.m). Duzeltme: esik
% HER SATIRIN KENDI b(i) buyuklugune gore ayri ayri hesaplaniyor.
n = length(b);
relTol = 1e-12;

if abs(b(1)) < relTol*max(1,abs(b(1)))
  warning('thomasSolve: near-zero pivot at row 1 (|b(1)|=%.3e)', abs(b(1)));
  x = nan(n,1);
  return;
end

cp = zeros(n,1); dp = zeros(n,1);
cp(1) = c(1)/b(1);
dp(1) = d(1)/b(1);
for i = 2:n
  m = b(i) - a(i)*cp(i-1);
  pivotTol_i = relTol * max(1, abs(b(i)));
  if abs(m) < pivotTol_i
    warning('thomasSolve: near-zero pivot at row %d (|m|=%.3e < tol=%.3e, |b(i)|=%.3e)', i, abs(m), pivotTol_i, abs(b(i)));
    x = nan(n,1);
    return;
  end
  cp(i) = c(i)/m;
  dp(i) = (d(i) - a(i)*dp(i-1))/m;
end
x = zeros(n,1);
x(n) = dp(n);
for i = n-1:-1:1
  x(i) = dp(i) - cp(i)*x(i+1);
end
end
