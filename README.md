# GaAsPN QFP drift-diffusion solver: publication snapshot

Version 1.0.0 (2026-11-19). Author: Volkan Kızıloğlu.
Repository: https://github.com/volkankiziloglu-dev/GaAsPN-DD-solver-publication

MATLAB code of a one-dimensional quasi-Fermi-potential (QFP) drift-diffusion solver for
heterojunction solar cells and the scripts that produce the results, figures and tables of the paper
comparing GaAsPN (Si-lattice-matched dilute nitride) solar-cell architectures (homojunction HJ-0,
front heterojunction F-HJ, GaP rear surface field RBSF-GaP, rear heterojunction RHJ, GaP/GaAsPN/GaP
p-i-n) under AM1.5G and indoor (CIE LED/fluorescent) illumination, together with the uncalibrated
validation against the GaInP cells of Klitzke et al. (2025).

Solver: Poisson and electron/hole continuity in QFP form, Scharfetter-Gummel fluxes with
heterointerface band offsets, SRH (one or two centres), radiative and Auger recombination, interface
recombination, Robin rear contact, Gummel-alternating block Newton iteration with analytic Jacobians,
TMM/Beer-Lambert optical generation with a two-layer ARC; acceptance gates for AM1.5G and low-current
indoor voltage points; optional F_1/2 Fermi-Dirac statistics, thermionic emission, and a double-double
/ excess-carrier refinement (`src/experimental_lowlight/`).

## Layout

| Folder / file | Content |
|---|---|
| `startup.m` | adds `src/`, `paper_scripts/`, `examples/` to the path and checks the solver location |
| `src/` | solver: 73 MATLAB files of dd_solver v1.1.2 selected by dependency analysis plus `resolveDataFile.m` (new), and 5 files in `src/experimental_lowlight/` (extended-precision refinement of the Supplementary Material); `src/data/` = redistributable spectra (CIE, ASTM G173) |
| `paper_scripts/` | the scripts that produce the paper results (flat folder; grouping by paper item in `docs/MANIFEST.md`); `paper_scripts/experimental_lowlight/` = drivers of the refinement |
| `external_validation/` | Sesame and SCAPS-1D cross-check drivers (Python) and their output summaries |
| `examples/` | two short examples (AM1.5G and 1000 lx J-V of RBSF-GaP) |
| `expected_outputs/results/` | summary CSVs of the paper runs (inputs of the figure/table scripts), 5.5 MB |
| `data/` | nominal parameter table and parameter-source matrix; `data/external/` for third-party files you obtain yourself |
| `docs/MANIFEST.md` | paper figure/table/section -> script -> inputs/outputs -> runtime |
| `docs/data_sources.md` | sources, citations and licences of all input data, and how to obtain the files that are not redistributed |

Generated output goes to `results/`, `figures/main_results/`, `latex/tables/` and `reports/`
(ignored by git).

## Requirements

- MATLAB R2024a (other releases not tested).
- Toolboxes (from `matlab.codetools.requiredFilesAndProducts`): none for the solver, the examples and
  almost all paper scripts. Statistics and Machine Learning Toolbox for
  `paper_scripts/analyze_uncertainty_sampling_ns.m`; Symbolic Math Toolbox for the variable-precision
  check of `paper_scripts/experimental_lowlight/` (methods `'prod+vpa'`, `'prod+vpaeval'`).
- Python 3 with numpy, pandas, scipy, matplotlib for the Python post-processing scripts; Sesame 2.0
  and SCAPS-1D 3.3.12 only for the external cross-checks (`docs/data_sources.md`).

## Third-party input data

The GaAsPN, GaP and GaPN optical tables and the digitized Klitzke spectrum/J-V curves are **not**
redistributed (unclear terms). Before running any device calculation, obtain them as described in
`docs/data_sources.md` and put them into `data/external/` or point the environment variable
`GAASPN_EXTERNAL_DATA` to a folder that contains them. The examples need only `nk_3N.mat`,
`GaP_Aspnes1983_nk.csv` and `GaPN_Shockhovets2018_digitized.csv`.

## Quick start

```matlab
% start MATLAB in the package folder (startup.m runs automatically), or:
cd path/to/GaAsPN-DD-solver-publication; startup
run examples/example_am15g_rbsf_gap.m     % RBSF-GaP, AM1.5G:   PCE 22.36 %  (~25 s)
run examples/example_indoor_1000lx.m      % RBSF-GaP, 6500 K, 1000 lx: PCE 45.45 %, PASS-A  (~12 s)
```

Smoke test (MATLAB R2024a, Windows 11, fresh session started in the package folder):

| Example | Printed metrics | Time |
|---|---|---|
| `example_am15g_rbsf_gap` | Jsc 16.815 mA/cm2, Voc 1467.1 mV, FF 90.66 %, PCE 22.36 % (79/79 points accepted) | 23 s (sweep 6 s) |
| `example_indoor_1000lx` | Pin 3.275 W/m2, class PASS-A, Jsc 129.22 uA/cm2, Voc 1324.9 mV, FF 86.95 %, PCE 45.45 % | 11 s (sweep 7 s) |

Both agree with the values of the paper (Table "AM1.5G performance"; indoor production run).

## Reproducing figures and tables

```matlab
restore_expected_outputs        % copy expected_outputs/results -> results (inputs of the figure/table scripts)
generate_heatmap_figures        % Figs. 7 and 9
generate_outdoor_jv_1col        % Fig. 5
generate_spatial_comparison_figure
generate_lux_pce_figure
```
```bash
python paper_scripts/plot_design_map.py          # Fig. 10
python paper_scripts/make_supp_tables_checks7.py # Supplementary tables -> latex/tables/
```

The device runs that create these summaries are listed per paper item in `docs/MANIFEST.md`
(e.g. `run_paper_indoor_production`, `run_paper_design_map`, `run_klitzke_table_v112`); most take
minutes to hours. Many runs read `results/paper_stageA/stageA_comparison.mat` (ARC thicknesses),
which is restored by `restore_expected_outputs` or recreated by `run_paper_stageA_comparison`.
Several scripts run the AM1.5G rows when the environment variable `AM_UNIFIED=1` is set and the
indoor rows otherwise; `merge_unified.py` combines them.

## Notes

- The paper scripts are copies of the project scripts; only paths were adapted (see `CHANGELOG.md`).
  Some comments are in Turkish.
- Tables 3 and 5 of the main text and a few Supplementary tables were assembled by hand from the CSV
  outputs named in `docs/MANIFEST.md`.

## Licence and citation

Code: BSD 3-Clause (`LICENSE`). Author-generated data and documentation: CC BY 4.0
(`LICENSE-data.md`); third-party files keep their own terms (CIE files: CC BY-SA 4.0).
Please cite the software as given in `CITATION.cff` (the archived release has a Zenodo DOI) and the paper.
