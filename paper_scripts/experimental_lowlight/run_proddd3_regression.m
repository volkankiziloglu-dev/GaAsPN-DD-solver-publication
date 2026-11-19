function run_proddd3_regression(part, stepTol)
% stepTol (optional): polishStateHPv3 flag threshold [V]; if given, outputs get the '_tol<value>' tag.
%RUN_PRODDD3_REGRESSION EXPERIMENTAL: 'prod+dd3' = prod+dd with polishStateHPv3 (new convergence-flag definition).
% part: 'm1' (100/200, 1000/2000/5000 lx, AM1.5G; mesh x1) | 'm2' | 'm4' (100, 1000 lx, AM1.5G). Production states are available under states/.
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
