# Data sources

## 1. Redistributed with the package (`src/data/`)

| File | Content | Source and citation | Terms |
|---|---|---|---|
| `CIE_illum_LEDs_1nm.csv` | Relative SPDs of the CIE LED illuminants (LED-B1 ... LED-V2), 380-780 nm, 1 nm; the paper uses LED-B1 (2700 K), LED-B3 (4000 K), LED-B5 (6500 K) | CIE, *Relative spectral power distributions of illuminants representing typical LED lamps, 1 nm spacing*, CIE 015:2018 dataset, DOI 10.25039/CIE.DS.dhcw57sd | CC BY-SA 4.0 (CIE datasets); unchanged |
| `CIE_illum_FLs_1nm.csv` | Relative SPDs of the CIE fluorescent illuminants FL1-FL12, FL3.1-FL3.15, 1 nm; the paper uses FL11 ("4000F") | CIE, *Relative spectral power distributions of illuminants representing typical fluorescent lamps, 1 nm wavelength steps*, CIE 015:2018 dataset, DOI 10.25039/CIE.DS.54hy6srn | CC BY-SA 4.0 (CIE datasets); unchanged |
| `CIE_sle_photopic.csv` | CIE photopic luminous efficiency V(λ), 360-830 nm, 1 nm | CIE spectral luminous efficiency for photopic vision, CIE datasets portal (https://cie.co.at/data-tables) | CC BY-SA 4.0 (CIE datasets); unchanged |
| `AM15.xls` | ASTM G173-03 reference spectra (extraterrestrial, global tilt AM1.5G, direct+circumsolar), 280-4000 nm | ASTM G173-03, *Standard Tables for Reference Solar Spectral Irradiances*; Excel file as distributed by NREL (https://www.nrel.gov/grid/solar-resource/spectra-am1.5.html) | Distributed freely by NREL; the standard is © ASTM International. Cite ASTM G173-03. |
| `AM15D.csv` | Direct+circumsolar (AM1.5D) column of the same table, two columns (nm, W m-2 nm-1), no header | as above (taken from the pvlib-python bundled copy of the NREL file) | as above |

The CC BY-SA 4.0 terms apply to the CIE files only (attribution: CIE; share-alike for
adaptations of these files); they do not extend to the code.

## 2. Not redistributed: obtain separately

These files are read through `src/resolveDataFile.m`. Put them in `data/external/` (or in the folder
named by the environment variable `GAASPN_EXTERNAL_DATA`) under the file names given here.
They were excluded because their redistribution terms are unclear (tables from journal articles, data
digitized from published figures, or data received from authors).

### 2.1 Needed for every device calculation (including the two examples)

| File (name expected by the code) | Content and format | Source / citation | How to obtain or recreate |
|---|---|---|---|
| `nk_3N.mat` (also `nk_25N.mat`, `nk_35N.mat`, `nk_4N.mat`, `nk_45N.mat`, `nk_5N.mat`) | Variable `nk_data` = [λ (µm), n, k] of Si-lattice-matched GaN_xP_yAs_(1-x-y) with x = 3 % (2.5-5 %) N. The paper uses 3 % N (nominal) and 2.5 % N (`run_paper_n25_comparison.m`). | Y. Zou, S. M. Goodnick, *Calculation of optical response functions of dilute-N GaPAsN lattice-matched to Si*, J. Appl. Phys. 127 (2020) 075703, DOI 10.1063/1.5140482 | Digitize the n(E), k(E) curves of the article or request the tables from the authors; save as `.mat` with variable `nk_data` (λ in µm, ascending). |
| `GaP_Aspnes1983_nk.csv` | Header `lambda_um,n,k`; GaP optical constants 0.21-0.83 µm (1.5-6.0 eV) | D. E. Aspnes, A. A. Studna, Phys. Rev. B 27 (1983) 985-1009, DOI 10.1103/PhysRevB.27.985 | Transcribe the GaP table of the article (also tabulated in the Handbook of Optical Constants of Solids). |
| `GaPN_Shockhovets2018_digitized.csv` | Header `E_eV,n,alpha_cm1`; GaP(N) lattice-matched to Si, 0.6-3.5 eV (used for the GaNP layer optics and for B_rad of GaNP in `src/paperRadiativeCoefficients.m`) | S. Shokhovets, O. Supplie, C. Koppka, S. Krischok, T. Hannappel, Phys. Rev. B 98 (2018) 075205, DOI 10.1103/PhysRevB.98.075205 | Digitize the n and α curves of the article. |

With these three files (plus the shipped CIE/ASTM files) both examples run; this was checked in the
smoke test.

### 2.2 Needed only for the Klitzke validation, the CW-LED rows and (AlInP) the Sukrittanon cells

| File | Content and format | Source / citation | How to obtain or recreate |
|---|---|---|---|
| `Klitzke2025_CWLED_SPD_digitized.csv` | Header `lambda_nm,relative_irradiance,abs_uncertainty`; cold-white LED spectrum, 380-780 nm, 1 nm | M. Klitzke et al., *Optimization of GaInP absorber design for indoor photovoltaic conversion efficiency above 40 %*, Appl. Phys. Lett. 127 (2025) 023301 (spectrum figure) | Digitize the measured LED spectrum of the article (shape only; the code scales it to 0.31 / 3.09 W m-2). |
| `experimental/klitzke2025/klitzke2025_fig4_iv_curves.csv` | Header `device,lux,V_V,J_uAcm2`; digitized J-V curves of the Klitzke HJ and RHJ devices at 100-1000 lx (read by `generate_klitzke_jv_1col.m`, `run_klitzke_rhj_excess_v112.m`, `fit_klitzke_hj_branch_v112.m`) | Klitzke et al. (2025), Fig. 4 | Digitize Fig. 4 of the article (e.g. with WebPlotDigitizer). |
| `AlInP_Ochoa2018_Cheong2015_optical.csv` | Header `wavelength_nm,n,alpha_cm1`; AlInP window optics | E. Ochoa-Martínez et al., Sol. Energy Mater. Sol. Cells 174 (2018) 388-396, DOI 10.1016/j.solmat.2017.09.028; J. S. Cheong et al., J. Phys. D 48 (2015) 405101, DOI 10.1088/0022-3727/48/40/405101 | Compile n from Ochoa-Martínez et al. and α near the edge from Cheong et al. |

### 2.3 Not used by the package

`schubert_data.mat` (GaInP optical constants), `AlInP_Cheong2015_digitized.csv` and the files of
`data/experimental/geisz2002/` and `data/experimental/sukrittanon2015/` of the project are not read
by any script of this package and are not needed.

### 2.4 Published values embedded in scripts

Several validation scripts contain a few published device metrics as numbers (measured J_sc, V_oc,
FF, PCE of the da Silva, Sukrittanon, Geisz MD421 and Klitzke cells; eye-digitized dark-J-V points of
the Sukrittanon thesis in `analyze_sukrittanon_*.m`), with the citation in the script header.
These are short factual values quoted for comparison.

### 2.5 Software not distributed

- **SCAPS-1D 3.3.12** (M. Burgelman et al., University of Gent), free for research on request from
  the authors. `external_validation/scaps_*.py` need the executable (`SCAPS_EXE`), a work folder
  (`SCAPS_WORK`) and a SCAPS definition file used as the text template (`SCAPS_TEMPLATE`; any
  three-layer `.def` file with one SRH defect in the last layer, as described in `scaps_build.py`).
- **Sesame 2.0** (B. Gaury et al., Sol. Energy Mater. Sol. Cells 198 (2019) 53), Python package
  (`pip install sesame`), BSD licence; needed by `external_validation/sesame_*.py`.

## 3. Author-generated data (`data/`, `expected_outputs/`, cross-check outputs)

- `data/paper_nominal_parameters.csv`: nominal parameters of the paper with status and basis.
- `data/material_parameter_source_matrix.csv`: literature sources and evidence tiers of the parameters.
- `expected_outputs/results/`: summary CSVs of the simulations (and `paper_stageA/stageA_comparison.mat`,
  the ARC thicknesses used by all runs). Some Klitzke-related files contain model-minus-measurement
  quantities derived from the digitized Klitzke curves
  (`paper_klitzke_rhj_excess_v112/rhj_excess_current_mesh4.csv`, used by the leakage analysis).
- `external_validation/scaps_output/`, `external_validation/sesame_output/`: cross-check J-V curves
  and summaries.

Licence: CC BY 4.0 (`LICENSE-data.md`).
