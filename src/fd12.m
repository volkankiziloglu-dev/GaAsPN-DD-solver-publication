function [F, dlnF, lnF] = fd12(eta)
%FD12 Normalized Fermi-Dirac integral of order 1/2, F(eta) = (2/sqrt(pi)) int_0^inf sqrt(x)/(1+exp(x-eta)) dx, so that n = Nc*F(eta)
% (2026-11-07, round-4 review R1-3: replaces the Blakemore approximation, which saturates at n/Nc = 1/gamma = 3.7 and fails in the
% accumulation layer at the n+-GaP/GaAsPN interface, eta = 6.4).
% Analytic approximation of Aymerich-Humet, Serra-Mestres and Millan (Solid-State Electron. 24 (1981) 981; also Blakemore 1982):
%   F = 1/(exp(-eta) + xi),  xi = 3*sqrt(pi/2) / u^(3/2),  u = eta + b + (|eta-b|^c + a^c)^(1/c),  with (j = 1/2)
%   a = sqrt(1 + 15/4 (j+1) + (j+1)^2/40) = 2.585, b = 1.8 + 0.61 j = 2.105, c = 2 + (2 - sqrt 2) 2^-j = 2.414 (a^c = 9.9),
% relative error below about 0.5 % for all eta, exact asymptotes exp(eta) (eta -> -inf) and 4/(3 sqrt(pi)) eta^(3/2) (eta -> +inf),
% smooth (C^1) and strictly increasing.
% Outputs: F, dlnF = d(ln F)/d(eta) (= 1/g, g the degeneracy factor of the generalized Einstein relation), lnF (stable).
j = 0.5; a = sqrt(1 + 15/4*(j+1) + (j+1)^2/40); b = 1.8 + 0.61*j; c = 2 + (2 - sqrt(2))*2^(-j); K = 3*sqrt(pi/2);
eta = max(min(eta, 700), -700);
s = eta - b; as = abs(s);
w = (as.^c + a^c);
u = eta + b + w.^(1/c);
xi = K * u.^(-1.5);
em = exp(-eta);
D = em + xi;
F = 1 ./ D;
lnF = -log(D);
du = 1 + w.^(1/c - 1) .* as.^(c - 1) .* sign(s);
dxi = -1.5 * K * u.^(-2.5) .* du;
dlnF = (em - dxi) ./ D;
end
