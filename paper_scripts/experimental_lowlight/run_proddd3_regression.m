function run_proddd3_regression(part, stepTol)
% stepTol (istege bagli): polishStateHPv3 bayrak esigi [V]; verilirse ciktilar '_tol<deger>' etiketini alir.
%RUN_PRODDD3_REGRESSION DENEYSEL: 'prod+dd3' = prod+dd, polishStateHPv3 (yeni yakinsama bayragi tanimi) ile.
% part: 'm1' (100/200, 1000/2000/5000 lx, AM1.5G; ag x1) | 'm2' | 'm4' (100, 1000 lx, AM1.5G). Uretim durumlari states/ altinda hazir.
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root,'paper_scripts','experimental_lowlight'));
base = struct(); if nargin >= 2, base.stepTol = stepTol; end
set_ = @(s) mergeOpts(base, s);
codes = {'HJ-0','F-HJ','RBSF-GaP','RHJ'};
switch part
  case 'm1'
    run_lowlight_experiment('prod+dd3', codes, [100 200], base);
    run_lowlight_experiment('prod+dd3', codes, [1000 2000 5000], base);
    run_lowlight_experiment('prod+dd3', codes, NaN, set_(struct('spectrum','AM15G')));
  otherwise
    ms = str2double(part(2:end));
    run_lowlight_experiment('prod+dd3', codes, 100, set_(struct('meshScale',ms)));
    run_lowlight_experiment('prod+dd3', codes, 1000, set_(struct('meshScale',ms)));
    run_lowlight_experiment('prod+dd3', codes, NaN, set_(struct('spectrum','AM15G','meshScale',ms)));
end
end

function o = mergeOpts(a, b)
o = a; for fn = fieldnames(b)', o.(fn{1}) = b.(fn{1}); end
end
