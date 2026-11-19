"""Analysis of results/paper_design_map_merged/design_map_*.csv: crossover ratio r*(tau_p) at which the RHJ equals RBSF-GaP, and the (tau_p, r) sign map.
r* is found by linear interpolation of PCE_RHJ - PCE_RBSF versus log10(r) between the grid ratios (RHJ increases monotonically with r);
r* < 1 means RHJ >= RBSF-GaP already at equal lifetimes, r* > 60 means no crossover in the grid. Also the HJ-0 comparison.
Output: results/paper_design_map_merged/design_map_summary.csv, rstar.csv; prints LaTeX-ready rows.
"""
import glob, os
import numpy as np, pandas as pd

root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
d = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(os.path.join(root, "results", "paper_design_map_merged", "design_map_[abcd].csv")))])
d.to_csv(os.path.join(root, "results", "paper_design_map_merged", "design_map_summary.csv"), index=False)
rows = []
for (cond, sb, tp), g in d.groupby(["Condition", "Sback", "tauP_ns"]):
    rb = g[g.Architecture == "RBSF-GaP"].iloc[0]; hj = g[g.Architecture == "HJ-0"].iloc[0]
    rh = g[g.Architecture == "RHJ"].sort_values("ratio")
    r = rh.ratio.values; diff = (rh.PCE_pct - rb.PCE_pct).values
    ok = np.isfinite(diff)
    r, diff = r[ok], diff[ok]
    if len(r) < 2: continue
    lr = np.log10(r)
    if diff[0] >= 0: rstar = 1.0 if diff[0] == 0 else np.nan   # RHJ ahead already at r = 1
    elif diff[-1] < 0: rstar = np.inf
    else:
        k = np.where(diff >= 0)[0][0]
        rstar = 10 ** (lr[k - 1] + (0 - diff[k - 1]) * (lr[k] - lr[k - 1]) / (diff[k] - diff[k - 1]))
    # RHJ at r=1 versus HJ-0, and the maximum possible RHJ gain over RBSF within the grid
    rows.append((cond, sb, tp, rb.PCE_pct, hj.PCE_pct, diff[0], diff[-1], rstar,
                 (rh.PCE_pct.values[0] - hj.PCE_pct), rb.Class, "/".join(sorted(set(rh.Class.astype(str))))))
R = pd.DataFrame(rows, columns=["Condition", "Sback", "tauP_ns", "PCE_RBSF", "PCE_HJ0", "RHJ_minus_RBSF_r1", "RHJ_minus_RBSF_r60", "r_star", "RHJ_r1_minus_HJ0", "Class_RBSF", "Classes_RHJ"])
R.to_csv(os.path.join(root, "results", "paper_design_map_merged", "rstar.csv"), index=False, float_format="%.3f")
pd.set_option("display.width", 220)
print(R.to_string(index=False, float_format=lambda v: "%.2f" % v))
