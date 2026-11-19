function R = tmmReflectance(lambda_nm, n_list, d_list_nm)
% Genel N-katmanli ince-film aktarim matrisi (transfer matrix method,
% Macleod "Thin-Film Optical Filters" formulasyonu) yansima hesaplayicisi,
% NORMAL GELIS (0 derece) icin. Herhangi bir ARC tasarimina uygulanabilir
% (Gorev #10), sadece GaP/RearHJT'ye ozgu degil.
%
% lambda_nm: dalga boyu [nm] (skaler)
% n_list: [n0 n1 n2 ... nK ns] -- kompleks kirilma indisleri, n0=gelis
%   ortami (genelde hava, 1.0), ns=alttas (yari-sonsuz), n1..nK ARA
%   katmanlar (ince film). Sogurma varsa n=n_gercek - i*k seklinde (fiziksel
%   isaret kurali: sonme, TMM formulunu degistirmez).
% d_list_nm: [d1 d2 ... dK] -- ARA katmanlarin kalinliklari [nm] (n_list'ten
%   2 eksik olmali: n0 ve ns yari-sonsuz oldugundan kalinliksiz).
%
% R: yansima orani (0-1), |r|^2.
n0 = n_list(1);
ns = n_list(end);
nFilms = n_list(2:end-1);
K = numel(nFilms);
if numel(d_list_nm) ~= K
  error('tmmReflectance: d_list_nm uzunlugu n_list''teki ARA katman sayisiyla (K=%d) eslesmiyor', K);
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
