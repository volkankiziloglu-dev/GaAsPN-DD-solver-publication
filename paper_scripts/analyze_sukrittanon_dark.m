function analyze_sukrittanon_dark()
%ANALYZE_SUKRITTANON_DARK Compare model dark J-V (tau_i = 5 ns, no leakage)
% with Sukrittanon thesis Fig. 4.7(a) (digitised by eye, ln J [mA/cm2], ~+-0.3).
root = fileparts(fileparts(mfilename('fullpath')));
rd = fullfile(root,'results','sukrittanon_ganp');
Rs = 22; % ohm cm2, from the light-J-V fit
% [V, ln(J/mA cm-2)]
ex.('GaNP2000') = [0.90 -7.3; 1.00 -6.0; 1.05 -5.3; 1.10 -4.4; 1.15 -3.2; 1.20 -2.0; 1.25 -1.2; 1.30 -0.55; 1.35 -0.3];
ex.('GaNP800')  = [0.60 -6.45; 0.80 -5.0; 0.90 -3.7; 1.00 -2.6; 1.10 -1.2; 1.20 -0.23; 1.30 1.35; 1.34 1.95];
ex.('GaP800')   = [0.60 -8.3; 0.70 -7.0; 0.80 -5.9; 0.90 -4.9; 1.00 -3.0; 1.10 -1.5; 1.20 0.3; 1.29 2.0];
names = {'GaNP2000','GaNP-2000'; 'GaNP800','GaNP-800'; 'GaP800','GaP-800'};
JscNoAR = [7.80 5.52 2.69]; % mA/cm2, measured, no AR
rows = [];
for i = 1:3
  m = readtable(fullfile(rd, sprintf('dark_%s_tau5e-09.csv', names{i,2})));
  e = ex.(names{i,1});
  Jm = interp1(m.V, log(m.Jdark_Acm2*1e3), e(:,1), 'linear', NaN);   % ln mA/cm2
  fprintf('%s: V, ln J exp, ln J model, exp/model\n', names{i,2});
  for k = 1:size(e,1)
    fprintf('  %.2f  %6.2f  %6.2f  x%.3g\n', e(k,1), e(k,2), Jm(k), exp(e(k,2)-Jm(k)));
  end
  % dark-derived Voc: Rs-corrected exp. dark curve at J = Jsc
  Jv = exp(e(:,2));                      % mA/cm2
  Vint = e(:,1) - Jv*1e-3*Rs;            % V
  [~,ia] = unique(e(:,2));
  Voc_dark = interp1(e(ia,2), Vint(ia), log(JscNoAR(i)), 'linear', NaN);
  fprintf('  dark-derived Voc (Rs=%d corrected, J_dark=Jsc=%.2f mA/cm2): %.3f V\n\n', Rs, JscNoAR(i), Voc_dark);
  rows = [rows; i Voc_dark]; %#ok<AGROW>
end
end
