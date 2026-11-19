function p = resolveDataFile(relName)
%RESOLVEDATAFILE Locate an input data file of the publication package.
%   p = resolveDataFile('CIE_illum_LEDs_1nm.csv')
%   p = resolveDataFile(fullfile('experimental','klitzke2025','klitzke2025_fig4_iv_curves.csv'))
%
% Search order (first hit wins):
%   1. <package>/src/data/<relName>        files redistributed with the package (CIE, ASTM G173)
%   2. <package>/data/external/<relName>   third-party files the user obtains separately
%   3. $GAASPN_EXTERNAL_DATA/<relName>     the same files kept in any other folder
%
% Third-party files that are not redistributed (optical constants, digitized spectra and J-V
% curves) are listed in docs/data_sources.md together with their sources.
srcDir = fileparts(mfilename('fullpath'));
pkgRoot = fileparts(srcDir);
cands = {fullfile(srcDir,'data',relName), fullfile(pkgRoot,'data','external',relName)};
ext = getenv('GAASPN_EXTERNAL_DATA');
if ~isempty(ext), cands{end+1} = fullfile(ext,relName); end
for k = 1:numel(cands)
  if isfile(cands{k}), p = cands{k}; return; end
end
error('resolveDataFile:NotFound', ['Input data file "%s" was not found.\n' ...
  'It is a third-party file that is not redistributed with this package. Obtain it as described in\n' ...
  'docs/data_sources.md and place it in data/external/ (same relative path) or in the folder named by\n' ...
  'the environment variable GAASPN_EXTERNAL_DATA.'], relName);
end
