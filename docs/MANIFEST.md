# Manifest: paper item -> script -> inputs/outputs

Paper: main text (`manuscript_en_kisa.tex`) and Supplementary Material (`supp_mat_en.tex`). Figure and table numbers follow the order of appearance in the LaTeX source.
All paths are relative to the package root. `results/` is created by the run scripts; the summary
files that the figure/table scripts need are shipped in `expected_outputs/results/` and are copied to
`results/` by `restore_expected_outputs` (MATLAB, in `paper_scripts/`).

Conventions

- **Run** = MATLAB device calculation (needs the third-party optical tables, `docs/data_sources.md`).
- **Post** = MATLAB or Python post-processing of saved CSVs (seconds).
- Runtime: one J-V curve at the nominal mesh (79 voltage points, 0.02 V) takes about 6-11 s on the
  test PC (MATLAB R2024a, single process); the 2x/4x meshes take roughly 2-4x and 4-10x longer.
  "n curves" gives the size of a run; times marked "measured" are from this package's smoke test or
  from timing columns of the original runs; all others are estimates (n curves x ~10 s).
- "Klitzke data" = digitized CW-LED spectrum and/or digitized J-V curves of Klitzke et al. (2025),
  not redistributed (see `docs/data_sources.md`).
- MATLAB scripts that are functions take the arguments shown in their header; with no arguments they
  run the full job used for the paper unless stated.

## Main text

| Paper item | Script(s) | Inputs | Outputs | Runtime |
|---|---|---|---|---|
| Fig. 1, model structure (TikZ) | none (hand-drawn `fig00_device_structure_en.tex`, not part of this package) | - | - | - |
| Fig. 2, indoor spectra (`fig_indoor_spectra`) | `generate_indoor_spectra_figure.m` (Post) | CIE files (shipped), Klitzke CW-LED SPD (Klitzke data) | `figures/main_results/fig_indoor_spectra.{pdf,png}` | s |
| Table 1, architectures (`table_architecture_hypotheses_kisa`) | none (descriptive table); device definitions in `src/buildPaperArchitecture.m` | - | - | - |
| Fig. 3, equilibrium band diagrams (`fig07_rear_band_alignment`) | `generate_rear_band_alignment_figure.m` (Run, equilibrium only) | optical tables not needed for equilibrium; `src/` | `fig07_rear_band_alignment.{pdf,png}`, `table_rear_band_offsets.csv` | 3.5 s (measured) |
| Table 2, material parameters (`paper_parameter_tables_en_kisa`) | none (hand-written); values in `src/buildPaperArchitecture.m`, `src/gapParams.m`, `data/paper_nominal_parameters.csv` | - | - | - |
| Fig. 4, Klitzke J-V model vs measurement (`fig_klitzke_jv_1col_stacked`) | `run_klitzke_hj_rhj_production.m` (Run; RHJ and HJ curves, `meshScale` 1 and 4), `run_klitzke_hj_bsf.m` (Run; HJ with BSF), then `generate_klitzke_jv_1col.m` (Post) | Klitzke data; `results/paper_klitzke_one_solver/*_mesh1_qfp.csv`, `results/paper_klitzke_hj_bsf/HJBSF_AlInP_*_curve.csv` | `figures/main_results/fig_klitzke_jv_1col_stacked.{pdf,png}` | runs: 4-6 curves each (~1-2 min at mesh 1); figure: s |
| Table 3, Klitzke model vs measurement (`table_klitzke_model_experiment_validation`) | `run_klitzke_hj_rhj_production.m`, `run_klitzke_table_v112.m` (AM1.5D rows, mesh 1/2/4), `run_klitzke_hj_bsf.m`, `fit_klitzke_hjbsf_terminal.m` (light-independent loss fit) | Klitzke data; ASTM G173 AM1.5D (`src/data/AM15D.csv`) | `results/paper_klitzke_one_solver/`, `results/paper_klitzke_table_v112/klitzke_table_v112_extra.csv`, `results/paper_klitzke_hj_bsf/terminal_fit.csv` | table assembled by hand from these CSVs (its former generator `generate_main_result_figures.m` writes outdated values and is not included) |
| Sec. "Code comparison" (Sesame, SCAPS-1D for Klitzke HJ/RHJ) | `export_klitzke_for_sesame_v112.m` (Run, export), `external_validation/sesame_klitzke_crosscheck.py`, `sesame_klitzke_summary_v112.py`, `scaps_klitzke_all.py` (+ `scaps_build_klitzke.py`, `scaps_run.py`) | Klitzke data; Sesame 2.0 (Python), SCAPS-1D 3.3.12 executable and a SCAPS definition template (not distributed) | `external_validation/sesame_output/klitzke_v112_crosscheck_summary.csv`, `external_validation/scaps_output/scaps_klitzke_v112_crosscheck.csv` (both shipped) | Sesame: minutes per device; SCAPS: s |
| Table 4, mesh convergence (`tab:mesh-convergence`) | Klitzke rows: `run_klitzke_table_v112.m` (HJ AM1.5D 2x; RHJ 100/1000 lx 4x), `run_klitzke_hj_rhj_production(4)` (HJ 100 lx 4x); HJ-0 row: `run_paper_mesh_richardson('HJ-0')` + `analyze_mesh_richardson.py` | Klitzke data (Klitzke rows) | `results/paper_klitzke_table_v112/*_mesh*.csv`, `results/paper_klitzke_one_solver/HJ_100lx_mesh4_qfp.csv`, `results/paper_mesh_richardson/richardson_summary.csv` | Richardson: 6 curves per architecture incl. 2x/4x, 168 s for RHJ (measured) |
| Sec. "Numerical convergence", voltage step, error budget | `run_paper_voltage_step_convergence.m` (AM1.5G 0.01 vs 0.02 V), `run_paper_indoor_current_1000lx.m` (1000 lx curves on 0.01 V grid), `test_unified_step.py`, `make_error_budget.py` | see SI rows below | `results/paper_voltage_step/`, `results/paper_mesh_richardson/unified_step_check.csv`, `latex/tables/supp_error_budget.tex` | runs: 12 + 16 curves |
| Fig. 5, AM1.5G J-V (`fig01_outdoor_jv_comparison_1col`) | `run_paper_stageA_comparison.m` (Run; ARC optimisation, writes `results/paper_stageA/stageA_comparison.mat` used by nearly all later runs), `run_paper_recombination_scenarios.m` (Run), `generate_outdoor_jv_1col.m` (Post) | AM1.5G (shipped) | `results/paper_recombination_scenarios/SRH+Rad+Auger_<arch>_jv.csv`, figure | scenarios: 28 curves (~5 min); figure 3 s (measured) |
| Table 5, AM1.5G performance (`table_outdoor_performance`) | `run_paper_recombination_scenarios.m` | - | `results/paper_recombination_scenarios/recombination_scenario_comparison.csv` (rows `SRH+Rad+Auger`) | table transcribed from the CSV (generator `generate_main_result_figures.m` not included) |
| Sec. 3.3 text: current balance at Voc | `run_recombination_decomposition.m` (Run) | - | `results/paper_recombination_decomposition/` | 3 curves |
| Fig. 6, spatial n(x), J_n(x) at MPP (`fig_spatial_comparison_MPP`) | `run_paper_spatial_profiles.m` (Run), `generate_spatial_comparison_figure.m` (Post) | - | `results/paper_spatial_profiles/<arch>_AM1.5G_MPP_{carriers,currents}.csv`, figure | 4 curves + profiles (~1 min); figure 2 s (measured) |
| Sec. 3.4 text: CW-LED results of the four designs | `run_paper_cwled_unified.m` (Run) | Klitzke CW-LED SPD | `results/paper_indoor_cwled_unified/` | 4 curves |
| Fig. 7, indoor PCE heat map (`fig03_indoor_pce_comparison`) | `run_paper_indoor_production.m` (Run; 7 architectures x 4 CIE spectra x 100/1000 lx), `make_indoor_unified_summary.py` (Post, unified extractor), `generate_heatmap_figures.m` (Post) | CIE files (shipped) | `results/paper_indoor_production/indoor_production_summary_unified.csv`, figure | 56 curves (~10 min); figure 12 s (measured, with Fig. 9) |
| Fig. 8, PCE vs illuminance (`fig_lux_pce`) | `run_paper_voc_lux.m` (Run, per architecture: `run_paper_voc_lux('RBSF-GaP')` etc.), `analyze_lux_series.py`, `generate_lux_pce_figure.m` (Post) | - | `results/paper_voc_lux/<arch>_tau<ns>_<lux>lx.csv`, `lux_series_pce.csv`, figure | 14 curves per architecture, 8-11 s each (measured) |
| Sec. "Recombination mechanisms" (radiative/Auger) | `run_paper_recombination_scenarios.m` | - | `recombination_scenario_comparison.csv` | see Fig. 5 |
| Sec. "Parameter robustness": permittivity / DOS | `run_paper_eps_dos_check.m` (Run; variants `'eps'`, `'dos'`) | - | `results/paper_eps_dos_check/` | 16 curves |
| Fig. 9, GaP chi x lifetime PCE maps (`fig04_gap_chi_lifetime_pce_sensitivity`) | `run_paper_gap_sensitivity_am.m` (Run), `generate_heatmap_figures.m` (Post) | - | `results/paper_gap_sensitivity_am/gap_sensitivity_summary.csv`, figure | not recorded (grid run) |
| Sec. "Absorber doping" | `run_paper_absorber_doping_ranking.m` (Run; AM1.5G rows with environment variable `AM_UNIFIED=1`), `merge_unified.py` | - | `results/paper_absorber_doping_unified[_am]/`, `results/paper_absorber_doping_merged/doping_ranking.csv` | 4 architectures x 3 dopings x 2 conditions per mode (plus unequal-doping rows) |
| Sec. "Paired parameter sampling" (main-text numbers; SI tables below) | `run_paper_uncertainty_sampling.m`, `run_paper_uncertainty_sampling_ns.m` (Run; set `AM_UNIFIED=1` for the AM1.5G rows), `merge_unified.py` | - | `results/paper_uncertainty_merged/uncertainty_samples.csv`, `results/paper_uncertainty_ns_merged/ns_samples_all.csv` | 480 curves per sampling (~1.5 h each) |
| Sec. "Further checks": Sesame / SCAPS for GaAsPN, band offset / thermionic emission | `export_for_sesame.m`, `export_for_sesame_fhj.m`, `export_for_scaps_prod.m` (Run, export), `external_validation/sesame_hetero_crosscheck.py`, `scaps_all.py`, `scaps_prod_all.py`; `run_paper_offset_thermionic.m` (Run) | Sesame 2.0; SCAPS-1D executable and template (not distributed) | `external_validation/sesame_output/*_sesame_JV.csv`, `external_validation/scaps_output/scaps_crosscheck.csv`, `scaps_prod_crosscheck.csv` (shipped); `results/paper_offset_thermionic/` | export: s; Sesame: minutes per device |
| Fig. 10, design map (`fig_design_map`) | `run_paper_design_map.m` (Run; indoor rows, then `AM_UNIFIED=1` for AM1.5G), `run_paper_design_map_mu.m` (Run, `muN` = 50/100/150), `merge_unified.py`, `analyze_design_map.py`, `analyze_rstar_band.py`, `analyze_design_map_mu.py`, `plot_design_map.py` (Post) | - | `results/paper_design_map_merged/{rstar,rstar_band}.csv`, `results/paper_design_map_mu/rstar_mu.csv`, figure | design map: several hours (not recorded); figure: s |
| Sec. "n-type absorber with a front junction" | `run_ntype_front_junction_check.m`, `run_ntype_fj_rstar.m` (Run), `analyze_ntype_fj_rstar.m` (Post) | - | `results/paper_ntype_front_junction/` | not recorded |
| Sec. "Nanosecond regime / calibrated regime" (MD421) | `run_geisz2002_staged_calibration.m` (Run, MD421 calibration), `run_paper_calibrated_regime.m` (Run, per architecture) | published MD421 metrics (embedded as numbers in the script) | `results/paper_calibrated_regime/` | finite-difference Jacobian solver, ~35 points per curve; minutes per architecture |
| Sec. "p-i-n" | `run_paper_pin.m` (Run; parts per header) | - | `results/paper_pin/pin_*.csv` | not recorded |
| Sec. "Unmodelled leakage" | `run_paper_leakage_robustness.m` (Post on production curves), `run_klitzke_rhj_excess_v112.m` (Run, RHJ excess current, mesh 4), `fit_klitzke_hj_branch_v112.m` (HJ branch fit), `analyze_leakage_v112.py` (Post) | Klitzke data | `results/paper_leakage_robustness/leakage_robustness.csv`, `results/paper_leakage_v112/leakage_archspecific_v112{,_wide}.csv` | post: s; excess-current run: 5 curves at mesh 4 |
| Sec. "Design guidance": shunt / perimeter numbers | `analyze_lux_series.py` (shunt variants), `analyze_perimeter_loss.py` | `results/paper_voc_lux/` | `results/paper_perimeter/` | s |
| Introduction: detailed-balance limits | `compute_sq_limits.m`, `compute_sq_cwled.m` (Post) | AM1.5G, CIE (shipped); CW-LED (Klitzke data) | `results/paper_sq_limits/` | s |

## Supplementary Material

| SI item | Script(s) | Inputs | Outputs | Runtime |
|---|---|---|---|---|
| Tables S (outdoor gates, indoor gates, curve classes) | none (descriptive); implemented in `src/qfpResultAccepted.m`, `src/indoorAcceptanceGate.m` | - | - | - |
| Sec. F-HJ curve classes at 1000 lx | `run_paper_indoor_production.m`; refinement: `paper_scripts/experimental_lowlight/run_lowlight_experiment.m` | - | `results/paper_indoor_production/` | see Fig. 7 |
| Table `supp_fd12_jv` (Fermi-Dirac F1/2 vs Boltzmann) | `run_paper_fd12_jv.m` (Run, per architecture and condition), `make_supp_table_fd12jv.py` (Post) | - | `results/paper_fd12_jv/*_summary.csv`, `fd12_jv_summary.csv`, `latex/tables/supp_fd12_jv.tex` | homotopy sweeps, not recorded |
| Sec. Bernoulli rounding | `run_paper_r4_analyses.m` (part `'BERN'`), `compare_bernoulli_fix.py` | - | `results/paper_r4_analyses/bernoulli_fix_comparison.csv` | 32 curves |
| Sec. Extended-precision refinement (double-double, excess-carrier flux) | `paper_scripts/experimental_lowlight/run_lowlight_experiment.m` (methods `'prod'`, `'prod+dd3'`, `'prod+vpa'`, ...), `run_proddd3_regression.m` (parts `'m1'`, `'m2'`, `'m4'`), `check_flux_forms.m`; solver files in `src/experimental_lowlight/` | production states are recomputed by method `'prod'`; `'prod+vpa'` needs the Symbolic Math Toolbox | `results/experimental_lowlight/summary_*.csv` (summaries shipped) | 48 curves plus refinement; `'prod+vpa'` is slow (hours) |
| Table `tab:supp-klitzke-devices` | none (device description); device builders inside the Klitzke run scripts | - | - | - |
| Table `tab:supp-hj-terminal` | `run_klitzke_hj_bsf.m`, `fit_klitzke_hjbsf_terminal.m` | Klitzke data; AM1.5D | `results/paper_klitzke_hj_bsf/` | 6 curves |
| Table `tab:supp-rhj-doping` | `run_klitzke_figs_v112('doping')` | Klitzke data | `results/paper_klitzke_figs_v112/rhj_doping_v112.csv` | not recorded |
| Table `tab:supp-rhj-excess` | `run_klitzke_rhj_excess_v112.m` | Klitzke data | `results/paper_klitzke_rhj_excess_v112/` | 5 curves at mesh 4 |
| Table `tab:supp-indoor-apriori` (`table_indoor_validation_status_en`) | none (hand-maintained summary of the Klitzke rows above) | - | - | - |
| Table `tab:supp-dasilva` | `run_dasilva_pin_qfp.m` | published da Silva metrics (numbers in the script); optical tables | `results/dasilva_pin_*` | 2 curves |
| Table `tab:supp-sukrittanon` | `run_sukrittanon_ganp_validation.m`, `run_sukrittanon_ganp_final.m`, `analyze_sukrittanon_dark.m`, `analyze_sukrittanon_leakage.m`, `analyze_sukrittanon_scaledSRH.m` | published Sukrittanon metrics and eye-digitized dark J-V points (embedded as numbers in the scripts; see `docs/data_sources.md`); GaNP optics (Shokhovets table) | `results/sukrittanon_ganp/` (selected CSVs shipped) | minutes |
| Table `tab:supp-voclux` | `run_paper_voc_lux.m`, `analyze_voc_lux.py` | - | printed table; `results/paper_voc_lux/` | see Fig. 8 |
| Tables `tab:supp-sesame`, `tab:supp-scaps`, `tab:supp-sesame-klitzke`; `supp_checks7_scaps_prod` | export scripts + `external_validation/*.py` (see main text rows); `make_supp_tables_checks7.py` writes `supp_checks7_scaps_prod.tex` | Sesame, SCAPS (not distributed) | `external_validation/*_output/` (shipped) | - |
| Table `tab:supp-mesh-gaaspn` | `run_paper_mesh_convergence.m` (HJ-U, RHJ; 2x mesh) | - | `results/paper_mesh_convergence/` | not recorded |
| Table `tab:supp-richardson`; mesh at 1 and 10 ns | `run_paper_mesh_richardson(code, scales, tauNs)`, `analyze_mesh_richardson.py` | - | `results/paper_mesh_richardson/richardson_summary.csv` | 168 s per architecture at 0.6 us (measured) |
| Sec. front contact, CV generation, heterointerface mesh | `run_paper_sfront.m`, `run_paper_mesh_cvgen.m`, `run_paper_interface_mesh.m` | - | `results/paper_sfront/`, `results/paper_mesh_cvgen/`, `results/paper_interface_mesh/` | not recorded |
| Table `supp_error_budget` | `make_error_budget.py` (Post) | `richardson_summary.csv`, `interface_mesh_*of6.csv`, `unified_step_check.csv`, `fd12_jv_summary.csv`, `bernoulli_fix_comparison.csv` | `results/paper_r4_analyses/error_budget.csv`, `latex/tables/supp_error_budget.tex` | s |
| Table `tab:supp-profiles`, Fig. S (RBSF-GaP profiles, `fig_spatial_RBSF_GaP_AM15G.png`) | `run_paper_spatial_profiles.m` (writes `figures/main_results/fig_spatial_RBSF_GaP_AM1.5G.{png,pdf}`; renamed to `..._AM15G.png` in the paper), Klitzke row: `run_klitzke_figs_v112('profiles')` | Klitzke data (Klitzke row) | `results/paper_spatial_profiles/`, `results/paper_klitzke_figs_v112/klitzke_rhj_1000lx_profiles_summary.csv` | ~1 min |
| Fig. S Klitzke RHJ profiles (`fig_klitzke_rhj_profiles`) | `run_klitzke_figs_v112('profiles')` (writes into `results/paper_klitzke_figs_v112/`) | Klitzke data | `fig_klitzke_rhj_profiles.{pdf,png}` | 1 curve |
| Table `tab:supp-luxseries` | `run_paper_voc_lux.m`, `analyze_lux_series.py` | - | `results/paper_voc_lux/lux_series_pce.csv` | s (post) |
| Tables `tab:supp-uncertainty-ranges/-pce/-order/-rho`, `tab:supp-ns-*`, `tab:supp-bootstrap`, `tab:supp-wtl` | `run_paper_uncertainty_sampling.m`, `run_paper_uncertainty_sampling_ns.m`, `merge_unified.py`, `analyze_uncertainty_sampling_ns.m` (Statistics and Machine Learning Toolbox), `bootstrap_spearman.py`, `analyze_win_tie_loss.py`, `make_sampling_tables_unified.py` (prints the table rows; argument `write` patches them into the LaTeX source, which is not part of this package) | - | `results/paper_uncertainty_merged/`, `results/paper_uncertainty_ns_merged/` | runs: ~1.5 h each; post: s |
| Tables `supp_checks7_prcc`, `supp_checks7_wtl_margin` | `analyze_sampling_r4_minor.py`, `make_supp_tables_checks7.py` | merged sampling CSVs | `r4_prcc.csv`, `r4_wtl_margins.csv`, `latex/tables/supp_checks7_*.tex` | s |
| Tables `tab:supp-edge`, `tab:supp-edge-eg`, `tab:supp-thick` | `run_paper_referee_checks2.m` (parts `'E'`, `'T'`), `run_paper_eg_optics_shift.m` (variants CB, VB), `merge_unified.py`, `analyze_eg_optics_shift.py`, `make_supp_tables_unified.py` (prints rows) | - | `results/paper_referee_checks_merged/` | not recorded |
| Table `tab:supp-n25` | `run_paper_n25_comparison.m` | 2.5 % N optical table (`nk_25N.mat`), CW-LED for one row | `results/paper_n25_comparison/` | not recorded |
| Table `tab:supp-offset` | `run_paper_offset_thermionic.m` | - | `results/paper_offset_thermionic/` | not recorded |
| Thermionic emission at high chi_GaP (Sec. r4B) | `run_paper_r4_analyses.m` (part per header), `analyze_r4_groupB.py` | - | `results/paper_r4_analyses/te_highchi.csv` | - |
| Rear-contact and lifetime-ratio checks; Table `tab:lifetime` | `run_paper_referee_checks.m` (part `'A'`), `run_paper_design_map.m`, `merge_unified.py`, `analyze_design_map.py` | - | `results/paper_referee_checks_merged/referee_checks_A.csv`, `results/paper_design_map_merged/` | see Fig. 10 |
| Table `supp_checks7_hju` (HJ-U control) | `run_paper_pin.m` (parts `'U'`, `'V'`), `make_supp_tables_checks7.py` | - | `results/paper_pin/pin_[UV]_*.csv` | - |
| Table `tab:supp-rstar` | `analyze_design_map.py` | design map | `results/paper_design_map_merged/rstar.csv` | s |
| Table `supp_dim_rstar` | `run_paper_design_map_mu.m`, `run_paper_r4_analyses.m` (thickness part), `analyze_dimensionless_criterion.py`, `make_supp_table_dim.py` | - | `results/paper_r4_analyses/dimensionless_rstar.csv`, `latex/tables/supp_dim_rstar.tex` | - |
| Table `supp_checks7_chi` | `run_paper_design_map_chi.m` (per chi_GaP), `analyze_design_map_chi.py`, `make_supp_tables_checks7.py` | - | `results/paper_design_map_chi/rstar_chi.csv` | - |
| Table `tab:supp-nfj` | `run_ntype_front_junction_check.m`, `run_ntype_fj_rstar.m`, `analyze_ntype_fj_rstar.m` | - | `results/paper_ntype_front_junction/` | 7 s (post, measured) |
| Table `supp_checks7_prior` | `analyze_prior_sensitivity.py`, `make_supp_tables_checks7.py` | `rstar.csv`, `rstar_mu.csv`, `ns_samples_all.csv` | `results/paper_prior_sensitivity/` | s |
| Table `tab:supp-calib` | `run_geisz2002_staged_calibration.m`, `run_paper_calibrated_regime.m` | - | `results/paper_calibrated_regime/` | see main text |
| Table `tab:supp-dasilva-tau` | `run_dasilva_tau_scan.m` | published J_sc values (in the script) | `results/paper_calibrated_regime/dasilva_tau_scan.csv` | 12 points |
| Tables `tab:supp-pin`, `tab:supp-pin-mu` | `run_paper_pin.m` | - | `results/paper_pin/` | - |
| Table `tab:supp-leakage-robustness` | `run_paper_leakage_robustness.m` (upper block), `analyze_leakage_v112.py` (architecture-specific block: rows `Set=production, HJbranch=old` of `leakage_archspecific_v112_wide.csv`) | Klitzke-derived excess current (`results/paper_klitzke_rhj_excess_v112/rhj_excess_current_mesh4.csv`) | table assembled from these CSVs | s |
| Table `supp_leak_archspecific_lux` | `analyze_leakage_v112.py` (rows `Set=voc_lux, HJbranch=old`, column `RBSF_minus_RHJ`) | as above | `results/paper_leakage_v112/leakage_archspecific_v112_wide.csv` | s; **note:** the shipped `make_supp_table_leak_lux.py` reads the earlier `leakage_archspecific_lux.csv` (from `analyze_leakage_archspecific_lux.py`, Klitzke model of the previous round) and reproduces the earlier table values; the published values equal the v112 file (author decision pending) |
| Table `supp_r4B_200lx` (200 lx design map and rear-contact controls) | `run_paper_r4_analyses.m` (part `'L200'`), `analyze_r4_groupB.py`, `make_supp_tables_r4B.py` | - | `results/paper_r4_analyses/rstar_200lx.csv`, `sback_200lx.csv` | - |
| Table `supp_r4B_perimeter` | `analyze_perimeter_loss.py`, `make_supp_tables_r4B.py` | `lux_series_pce.csv` and curves of `results/paper_voc_lux/` | `results/paper_perimeter/` | s |
| Table `supp_r4B_srh` (one SRH centre in both absorbers) | `run_paper_r4_analyses.m` (part `'SRH2'`), `analyze_r4_groupB.py`, `make_supp_tables_r4B.py` | - | `results/paper_r4_analyses/srh2_*.csv` | - |
| Table `tab:supp-sq` (detailed balance) | `compute_sq_limits.m`, `compute_sq_cwled.m` | AM1.5G, CIE (shipped), CW-LED (Klitzke data) | `results/paper_sq_limits/` | s |

## Verified in this package (smoke test, MATLAB R2024a, Python 3 / pandas / scipy / matplotlib)

- `examples/example_am15g_rbsf_gap.m` and `examples/example_indoor_1000lx.m` reproduce the paper values
  (22.36 % and 45.45 %).
- From `expected_outputs/`: `generate_heatmap_figures`, `generate_outdoor_jv_1col`,
  `generate_spatial_comparison_figure`, `generate_lux_pce_figure`, `generate_rear_band_alignment_figure`,
  `analyze_ntype_fj_rstar`, `analyze_uncertainty_sampling_ns`, `plot_design_map.py`,
  `make_supp_table_fd12jv.py`, `make_error_budget.py`, `make_supp_tables_checks7.py`,
  `make_supp_tables_r4B.py`, `make_supp_table_dim.py`, `make_sampling_tables_unified.py`,
  `make_supp_tables_unified.py`, `analyze_voc_lux.py` all run. The generated LaTeX tables agree with
  the paper tables except for caption wording edited by hand in the paper (and
  `supp_leak_archspecific_lux.tex`, see the note above).
- The long device runs were not repeated for this package.
