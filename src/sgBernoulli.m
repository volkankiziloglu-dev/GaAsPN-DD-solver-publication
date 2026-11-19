function y = sgBernoulli(x)
% B(x) = x / (exp(x)-1), robust near x=0 and for large |x|
y = zeros(size(x));
small = abs(x) < 1e-6;
big_pos = x > 80;
big_neg = x < -80;
mid = ~small & ~big_pos & ~big_neg;

y(small) = 1 - x(small)/2 + (x(small).^2)/12;
y(big_pos) = x(big_pos) .* exp(-x(big_pos));
y(big_neg) = -x(big_neg);
y(mid) = x(mid) ./ (exp(x(mid)) - 1);
end
