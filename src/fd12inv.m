function eta = fd12inv(x)
%FD12INV Inverse of fd12: eta such that F_{1/2}(eta) = x (x > 0), by Newton iteration on ln F from the Joyce-Dixon-like start
% eta0 = ln x (non-degenerate) or (3 sqrt(pi) x / 4)^(2/3) (degenerate). Used for equilibrium initial guesses and quasi-Fermi diagnostics.
x = max(x, 1e-300);
eta = log(x);
k = x > 1; eta(k) = max(log(x(k)), (3*sqrt(pi)*x(k)/4).^(2/3));
lx = log(x);
for it = 1:60
  [~, dl, lF] = fd12(eta);
  d = (lF - lx) ./ dl;
  eta = eta - d;
  if max(abs(d)) < 1e-13, break; end
end
end
