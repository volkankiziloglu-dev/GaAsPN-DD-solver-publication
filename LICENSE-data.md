# Licence of data and documentation

## Author-generated material: CC BY 4.0

Unless stated otherwise below, the data and documentation in this repository that were
created by the author are licensed under the Creative Commons Attribution 4.0 International
licence (CC BY 4.0, https://creativecommons.org/licenses/by/4.0/):

- `expected_outputs/` (summary results of the simulations reported in the paper),
- `external_validation/scaps_output/`, `external_validation/sesame_output/` (cross-check summaries),
- `data/paper_nominal_parameters.csv`, `data/material_parameter_source_matrix.csv`,
- `README.md`, `CHANGELOG.md` and the files in `docs/`.

Attribution: Volkan Kızıloğlu, *GaAsPN QFP drift-diffusion solver: publication snapshot*
(see `CITATION.cff`).

The source code is licensed separately under the BSD 3-Clause licence (`LICENSE`).

## Third-party files: their own terms

Files obtained from third parties keep their original terms; they are **not** covered by the
CC BY 4.0 or BSD licences above.

| File(s) | Source | Terms |
|---|---|---|
| `src/data/CIE_illum_LEDs_1nm.csv` | CIE, relative SPDs of illuminants representing typical LED lamps, 1 nm, CIE 015:2018 dataset, DOI 10.25039/CIE.DS.dhcw57sd | CC BY-SA 4.0 (CIE datasets); redistributed unchanged |
| `src/data/CIE_illum_FLs_1nm.csv` | CIE, relative SPDs of illuminants representing typical fluorescent lamps, 1 nm, CIE 015:2018 dataset, DOI 10.25039/CIE.DS.54hy6srn | CC BY-SA 4.0 (CIE datasets); redistributed unchanged |
| `src/data/CIE_sle_photopic.csv` | CIE spectral luminous efficiency for photopic vision, V(λ), CIE dataset (CIE datasets portal, https://cie.co.at/data-tables) | CC BY-SA 4.0 (CIE datasets); redistributed unchanged |
| `src/data/AM15.xls`, `src/data/AM15D.csv` | ASTM G173-03 reference spectra (AM1.5G global tilt, AM1.5D direct+circumsolar) as distributed by NREL (https://www.nrel.gov/grid/solar-resource/spectra-am1.5.html); `AM15D.csv` is the direct+circumsolar column of the same table (via the pvlib-python bundled copy) | Freely distributed by NREL; the standard itself is © ASTM International. Cite ASTM G173-03. See `docs/data_sources.md`. |

`expected_outputs/results/paper_klitzke_rhj_excess_v112/rhj_excess_current_mesh4.csv` is
author-generated (CC BY 4.0) but derived from J-V curves digitized from M. Klitzke et al.,
Appl. Phys. Lett. 127 (2025) 023301, Fig. 4 (measured minus simulated recombination current);
cite that article when using it.

Files that could not be redistributed (optical constants and digitized curves) are listed in
`docs/data_sources.md` with instructions for obtaining them.
