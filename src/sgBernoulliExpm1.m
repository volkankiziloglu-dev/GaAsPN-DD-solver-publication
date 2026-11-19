function y = sgBernoulliExpm1(x)
% B(x) = x / (exp(x)-1) with expm1 in the middle range (accurate to rounding for small |x|). The production sgBernoulli.m keeps
% x./(exp(x)-1), whose relative error ~1e-16/|x| leaves equilibrium segment currents of up to ~1e-3 A/cm2 in heavily doped,
% finely meshed layers; switching production to expm1 makes the unscaled Gummel line search fragile (round-5 review, 2026-11-08,
% results/paper_r4_analyses/bernoulli_fix_comparison.csv). Used by the F_1/2 residual (assembleResidualQFPFD12.m).
y = zeros(size(x));
small = abs(x) < 1e-6;
big_pos = x > 80;
big_neg = x < -80;
mid = ~small & ~big_pos & ~big_neg;

y(small) = 1 - x(small)/2 + (x(small).^2)/12;
y(big_pos) = x(big_pos) .* exp(-x(big_pos));
y(big_neg) = -x(big_neg);
y(mid) = x(mid) ./ expm1(x(mid));
end
