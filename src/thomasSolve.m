function x = thomasSolve(a, b, c, d)
% Solve tridiagonal system a(i)x(i-1)+b(i)x(i)+c(i)x(i+1)=d(i). a(1),c(n) unused.
%
% FIX (external review item 1.3): PREVIOUSLY there was NO protection for a near-zero/small pivot
% (m = b(i)-a(i)*cp(i-1)) -- it could produce Inf/NaN/large cancellation
% errors (the rationale of the warm-up ramp noted that in the past exactly
% this kind of ill-conditioning had occurred, see
% MATHEMATICAL_MODEL.md Section K). If the absolute pivot value falls below a threshold
% normalized to EACH ROW'S OWN b(i) scale (row-local),
% a WARNING is issued and a NaN vector is returned.
%
% IMPORTANT: the first version used a SINGLE, GLOBAL threshold (norm(b,inf)
% based) -- this gave FALSE ALARMS on the boundary-condition rows (b(1)=1, a deliberately
% normalized Dirichlet equation), because the diffusion coefficients of interior
% nodes (b can be of order 1e13-1e15,
% see K.12's Case B trace) dominated the global norm and made b(1)=1
% look "near zero" (caught IMMEDIATELY when the cross-check script was run
% -- debug_thomassolve_crosscheck.m). Fix: the threshold is computed
% separately according to EACH ROW'S OWN b(i) magnitude.
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
