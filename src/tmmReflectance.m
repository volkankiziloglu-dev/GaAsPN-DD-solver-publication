function R = tmmReflectance(lambda_nm, n_list, d_list_nm)
% General N-layer thin-film transfer matrix (transfer matrix method,
% Macleod "Thin-Film Optical Filters" formulation) reflectance calculator,
% for NORMAL INCIDENCE (0 degrees). Applicable to any ARC design
% (Task #10), not specific to GaP/RearHJT.
%
% lambda_nm: wavelength [nm] (scalar)
% n_list: [n0 n1 n2 ... nK ns] -- complex refractive indices, n0=incidence
%   medium (usually air, 1.0), ns=substrate (semi-infinite), n1..nK INTERMEDIATE
%   layers (thin films). With absorption, n=n_real - i*k (physical
%   sign convention: extinction, does not change the TMM formula).
% d_list_nm: [d1 d2 ... dK] -- thicknesses of the INTERMEDIATE layers [nm] (must be
%   2 shorter than n_list: n0 and ns are semi-infinite, hence have no thickness).
%
% R: reflectance (0-1), |r|^2.
n0 = n_list(1);
ns = n_list(end);
nFilms = n_list(2:end-1);
K = numel(nFilms);
if numel(d_list_nm) ~= K
  error('tmmReflectance: d_list_nm length does not match the number of INTERMEDIATE layers in n_list (K=%d)', K);
end
M = eye(2);
for j = 1:K
  nj = nFilms(j); dj = d_list_nm(j);
  delta = 2*pi*nj*dj/lambda_nm;
  Mj = [cos(delta), 1i*sin(delta)/nj; 1i*nj*sin(delta), cos(delta)];
  M = M*Mj;
end
BC = M*[1; ns];
B = BC(1); C = BC(2);
r = (n0*B - C)/(n0*B + C);
R = abs(r)^2;
end
