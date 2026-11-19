function assertCanonicalPath()
% External review item 4.7: MATLAB path precedence (startup.m adding the
% ENTIRE matlab/ tree via genpath, risk of COLLISION with SAME-NAMED but eps0-buggy
% old solver files in RearHJT/) relied on alphabetical ordering
% (GaAsNP_RHJ_final < RearHJT) -- works, but a fragile
% guarantee. This function, called at the start of every script that produces
% main results, EXPLICITLY verifies that the CORRECT (this folder's) solver
% files are actually used -- if an old RearHJT/ file is loaded by mistake,
% it stops immediately with a clear error instead of SILENTLY producing
% wrong results.
% FIX (external review item 1.1): startsWith(p,here) could THEORETICALLY
% match by mistake a different folder with the same PREFIX as "here" (e.g.
% GaAsNP_RHJ_final2/). The expected root including the folder BOUNDARY
% (ending with filesep) is used instead.
here = fileparts(mfilename('fullpath'));
expectedRoot = [here filesep];
fns = {'newtonPoisson','newtonPoissonML','gummelSolve','gummelSolveML', ...
       'sweepJV','sweepJVML','setupDevice','setupDeviceML','thomasSolve', ...
       'gummelSolveQFPAnalytic','gummelSolveQFP','assembleResidualQFP','interfaceRecombination', ...   % dd_solver v1.1.1: the QFP production path
       'runAcceptedJVSweep','indoorAcceptanceGate','pvMetrics','buildDeviceML', ...
       'sgBernoulli','qfpResultAccepted'};
% Publication package: legacy solver files (gummelSolve, sweepJV, setupDevice) are not part of the
% dependency closure and are therefore not shipped; only the functions present in src/ are checked.
fns = fns(cellfun(@(f) isfile(fullfile(here,[f '.m'])), fns));
for i = 1:numel(fns)
  p = which(fns{i});
  if isempty(p)
    error('assertCanonicalPath: %s not found at all (it may not have been added to the path)', fns{i});
  end
  if ~startsWith(p, expectedRoot)
    error(['assertCanonicalPath: %s is being loaded from an UNEXPECTED location:\n  %s\n' ...
      'Expected folder: %s\n' ...
      'This indicates that the old (eps0-buggy) RearHJT/ copy won path precedence -- ' ...
      'results may be SILENTLY wrong. Check the addpath order.'], ...
      fns{i}, p, here);
  end
end
end
