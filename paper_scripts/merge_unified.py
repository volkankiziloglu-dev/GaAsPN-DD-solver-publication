"""Merge AM1.5G rows and indoor rows, both with the unified extractor (all finite points), into the layouts used by the analysis scripts:
 results/paper_uncertainty_merged/uncertainty_samples.csv
 results/paper_uncertainty_ns_merged/ns_samples_all.csv
 results/paper_design_map_merged/design_map_[abcd].csv
 results/paper_referee_checks_merged/*.csv
 results/paper_absorber_doping_merged/doping_ranking.csv
AM1.5G rows: from the *_am directories (run with AM_UNIFIED=1; metrics from all finite points); if a *_am file is missing the earlier gate-accepted AM1.5G rows are used.
Indoor rows: from the *_unified directories."""
import os, glob
import pandas as pd
root = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); R = os.path.join(root, "results")
P = lambda *a: os.path.join(R, *a)


def read_all(files):
    return pd.concat([pd.read_csv(f) for f in files], ignore_index=True) if files else None


def merge(old_files, am_files, indoor_files, out_file, cond_col="Condition", sortcols=("Sample", "Architecture")):
    old = read_all(old_files); old_am = old[old[cond_col] == "AM1.5G"]
    am = read_all(am_files)
    am = am[am[cond_col] == "AM1.5G"] if am is not None else old_am
    ind = read_all(indoor_files); ind = ind[ind[cond_col] != "AM1.5G"]
    cols = list(old.columns)
    m = pd.concat([am[cols], ind[cols]], ignore_index=True)
    keys = [c for c in sortcols if c in cols] + [cond_col]
    m = m.sort_values(keys, kind="stable").reset_index(drop=True)
    os.makedirs(os.path.dirname(out_file), exist_ok=True); m.to_csv(out_file, index=False)
    print(os.path.relpath(out_file, R), "AM", len(am), "(old %d)" % len(old_am), "indoor", len(ind), "total", len(m))


merge([P("paper_uncertainty", "uncertainty_samples.csv")], sorted(glob.glob(P("paper_uncertainty_unified_am", "uncertainty_samples_*_*.csv"))),
      [P("paper_uncertainty_unified", "uncertainty_samples.csv")], P("paper_uncertainty_merged", "uncertainty_samples.csv"))
merge([P("paper_uncertainty_ns", "ns_samples_all.csv")], sorted(glob.glob(P("paper_uncertainty_ns_unified_am", "ns_samples_[0-9]*_[0-9]*.csv"))),
      sorted(glob.glob(P("paper_uncertainty_ns_unified", "ns_samples_[0-9]*_[0-9]*.csv"))), P("paper_uncertainty_ns_merged", "ns_samples_all.csv"))
for t in "abcd":
    fo = P("paper_design_map", f"design_map_{t}.csv"); fa = P("paper_design_map_unified_am", f"design_map_{t}.csv"); fn = P("paper_design_map_unified", f"design_map_{t}.csv")
    if os.path.exists(fn):
        merge([fo], [fa] if os.path.exists(fa) else [], [fn], P("paper_design_map_merged", f"design_map_{t}.csv"), sortcols=("tauP_ns", "Sback", "ratio", "Architecture"))
merge([P("paper_absorber_doping", "doping_ranking.csv")], [P("paper_absorber_doping_unified_am", "doping_ranking.csv")] if os.path.exists(P("paper_absorber_doping_unified_am", "doping_ranking.csv")) else [],
      [P("paper_absorber_doping_unified", "doping_ranking.csv")], P("paper_absorber_doping_merged", "doping_ranking.csv"), sortcols=("N_cm3", "Architecture"))

# ---- referee checks / Eg-optics shift
od = P("paper_referee_checks"); nd = P("paper_referee_checks_unified"); ad = P("paper_referee_checks_unified_am"); md = P("paper_referee_checks_merged")
os.makedirs(md, exist_ok=True)
for fn in ("referee_checks_A.csv", "referee_checks_B.csv", "referee_checks_C.csv", "referee_checks_E.csv", "referee_checks_T.csv", "eg_optics_shift_CB.csv", "eg_optics_shift_VB.csv"):
    if not os.path.exists(os.path.join(nd, fn)): continue
    o = pd.read_csv(os.path.join(od, fn)); n = pd.read_csv(os.path.join(nd, fn))
    a = pd.read_csv(os.path.join(ad, fn)) if os.path.exists(os.path.join(ad, fn)) else o
    a = a[a.Condition == "AM1.5G"]; n = n[n.Condition != "AM1.5G"]
    m = pd.concat([a[o.columns], n[o.columns]], ignore_index=True)
    m.to_csv(os.path.join(md, fn), index=False); print(fn, "AM", len(a), "(old %d)" % len(o[o.Condition == "AM1.5G"]), "indoor", len(n), "total", len(m))
