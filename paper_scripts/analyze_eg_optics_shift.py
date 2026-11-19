"""Summary of results/paper_referee_checks/eg_optics_shift_{CB,VB}.csv together with the earlier optical-only shift (referee_checks_E.csv): PCE of RBSF-GaP, ranking and margins versus shift."""
import os
import numpy as np, pandas as pd
root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
d = os.path.join(root, "results", "paper_referee_checks_merged")
cb = pd.read_csv(os.path.join(d, "eg_optics_shift_CB.csv")); vb = pd.read_csv(os.path.join(d, "eg_optics_shift_VB.csv"))
old = pd.read_csv(os.path.join(d, "referee_checks_E.csv")); old["variant"] = "optics only"; old["BradAbs"] = np.nan
allv = pd.concat([cb, vb, old[["Architecture", "shift_eV", "variant", "BradAbs", "Condition", "Jsc", "Voc_mV", "FF_pct", "PCE_pct", "Class"]]])
allv.to_csv(os.path.join(d, "eg_optics_shift_all.csv"), index=False, float_format="%.4f")
pd.set_option("display.width", 220)
rows = []
for (var, cond, sh), g in allv.groupby(["variant", "Condition", "shift_eV"]):
    p = g.set_index("Architecture").PCE_pct
    cls = g.set_index("Architecture").Class
    rows.append((var, cond, sh, p["HJ-0"], p["F-HJ"], p["RBSF-GaP"], p["RHJ"], p["RBSF-GaP"] - p["RHJ"], p["RBSF-GaP"] - p["HJ-0"],
                 " ".join(p.sort_values(ascending=False).index.str.replace("RBSF-GaP", "RB").tolist()), g.BradAbs.iloc[0]))
R = pd.DataFrame(rows, columns=["variant", "condition", "shift_eV", "HJ-0", "F-HJ", "RBSF", "RHJ", "RBSF-RHJ", "RBSF-HJ0", "order", "BradAbs"])
R.to_csv(os.path.join(d, "eg_optics_shift_summary.csv"), index=False, float_format="%.4f")
print(R.to_string(index=False, float_format=lambda v: "%.3g" % v if abs(v) < 1e-3 else "%.2f" % v))
