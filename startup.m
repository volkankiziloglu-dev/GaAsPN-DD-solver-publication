function startup()
%STARTUP Path setup of the GaAsPN QFP drift-diffusion publication package.
% Run from the package root (MATLAB runs it automatically when started in this folder).
% Adds src/, paper_scripts/ and examples/ and checks that the solver is loaded from src/.
% src/experimental_lowlight/ (extended-precision refinement) is added by the scripts that use it.
root = fileparts(mfilename('fullpath'));
addpath(fullfile(root,'examples'), '-begin');
addpath(fullfile(root,'paper_scripts'), '-begin');
addpath(fullfile(root,'src'), '-begin');

required = {'gummelSolveQFPAnalytic','buildPaperArchitecture','deviceGenAndTransmitML', ...
    'runAcceptedJVSweep','indoorAcceptanceGate','pvMetrics','resolveDataFile'};
expected = [fullfile(root,'src') filesep];
for k = 1:numel(required)
    p = which(required{k});
    if isempty(p)
        error('startup:MissingDependency', 'Required function not found: %s', required{k});
    end
    if ~startsWith(p, expected)
        error('startup:PathCollision', '%s is loaded from an unexpected location: %s', required{k}, p);
    end
end

ext = fullfile(root,'data','external');
if ~isfile(fullfile(ext,'nk_3N.mat')) && isempty(getenv('GAASPN_EXTERNAL_DATA'))
    fprintf(['Note: third-party optical data (e.g. nk_3N.mat) not found in data/external/ and GAASPN_EXTERNAL_DATA is not set.\n' ...
             '      Device calculations need them; see docs/data_sources.md.\n']);
end
fprintf('GaAsPN QFP drift-diffusion package ready: %s\n', root);
end
