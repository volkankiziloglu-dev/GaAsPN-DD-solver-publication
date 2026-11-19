# Changelog

## 1.0.0 (2026-11-19)

Publication snapshot of the solver and the paper scripts.

- Solver: dd_solver v1.1.2 (synced 2026-11-13), reduced to the dependency closure of the paper
  scripts (`matlab.codetools.requiredFilesAndProducts`), plus `src/experimental_lowlight/`
  (double-double / excess-carrier refinement of the Supplementary Material).
- Changes made for the package (copies only; numerical code unchanged):
  - `src/resolveDataFile.m` (new): locates input data in `src/data/`, `data/external/` or
    `$GAASPN_EXTERNAL_DATA`; the optical-data and spectrum loaders and the paper scripts use it
    instead of fixed paths, because several third-party files are not redistributed.
  - `src/assertCanonicalPath.m`: checks only the solver functions that are shipped.
  - Paper scripts moved from `cases/` to `paper_scripts/`; `addpath` calls adapted;
    the `tests/` path of the project removed (not needed).
  - SCAPS helper scripts: executable, work folder and template file taken from the environment
    variables `SCAPS_EXE`, `SCAPS_WORK`, `SCAPS_TEMPLATE` instead of user-specific paths.
  - Local absolute paths removed from comments and `SOLVER_VERSION.txt`.
- Added: `startup.m`, two short examples, `paper_scripts/restore_expected_outputs.m`,
  documentation (`docs/`), licences and citation metadata.
