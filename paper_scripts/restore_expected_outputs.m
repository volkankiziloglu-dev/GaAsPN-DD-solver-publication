function restore_expected_outputs(overwrite)
%RESTORE_EXPECTED_OUTPUTS Copy the summary CSVs shipped in expected_outputs/results/ to results/.
% The figure and table scripts read their inputs from results/<folder>/. Running this function lets them be executed
% without first repeating the (hours to days of) device calculations. Existing files in results/ are kept unless
% overwrite = true. results/paper_stageA/stageA_comparison.mat (ARC thicknesses) is needed by almost every run script.
if nargin < 1, overwrite = false; end
root = fileparts(fileparts(mfilename('fullpath')));
src = fullfile(root,'expected_outputs','results'); dst = fullfile(root,'results');
L = dir(fullfile(src,'**','*.*')); L = L(~[L.isdir]);
n = 0;
for k = 1:numel(L)
  rel = erase(fullfile(L(k).folder, L(k).name), [src filesep]);
  out = fullfile(dst, rel);
  if isfile(out) && ~overwrite, continue; end
  if ~isfolder(fileparts(out)), mkdir(fileparts(out)); end
  copyfile(fullfile(L(k).folder, L(k).name), out); n = n + 1;
end
fprintf('%d files copied from expected_outputs/results to results/\n', n);
end
