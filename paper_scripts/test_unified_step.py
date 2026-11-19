"""Effect of the extractor on the voltage-step sensitivity of the indoor 1000 lx PCE: the 0.02 V production curves and the 0.01 V current-code curves, each evaluated with the
old indoor extractor (grid maximum, as stored in the summaries) and with the unified extractor (MPP parabola). Output: results/paper_mesh_richardson/unified_step_check.csv"""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from unified_metrics import pv_metrics
root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
S = pd.read_csv(os.path.join(root, "results", "paper_indoor_production", "indoor_production_summary.csv"))
C = pd.read_csv(os.path.join(root, "results", "paper_indoor_current_1000lx", "current_1000lx_summary.csv"))
out = []
for a in ["HJ-0", "F-HJ", "RBSF-GaP", "RHJ"]:
    for sp in ["2700K", "4000K", "6500K", "4000F"]:
        r = S[(S.Architecture == a) & (S.Spectrum == sp) & (S.Lux == 1000)].iloc[0]
        c = pd.read_csv(os.path.join(root, "results", "paper_indoor_production", f"{a}_{sp}_1000lux_jv.csv")); c = c[np.isfinite(c.J_Acm2)]
        m = pv_metrics(c.V_V, c.J_Acm2, r.Pin_Wm2)
        q = C[(C.Architecture == a) & (C.Spectrum == sp)].iloc[0]
        c2 = pd.read_csv(os.path.join(root, "results", "paper_indoor_current_1000lx", f"{a}_{sp}_1000lx.csv")); c2 = c2[np.isfinite(c2.J_Acm2)]
        m2 = pv_metrics(c2.V_V, c2.J_Acm2, q.Pin_Wm2)
        out.append((a, sp, r.PCE_pct, m["PCE"], q.PCE_pct, m2["PCE"]))
d = pd.DataFrame(out, columns=["Architecture", "Spectrum", "old_0.02V", "unified_0.02V", "old_0.01V", "unified_0.01V"])
d["old_step_diff"] = d["old_0.01V"] - d["old_0.02V"]; d["unified_step_diff"] = d["unified_0.01V"] - d["unified_0.02V"]
d.to_csv(os.path.join(root, "results", "paper_mesh_richardson", "unified_step_check.csv"), index=False, float_format="%.4f")
pd.set_option("display.width", 200)
print(d.to_string(index=False, float_format=lambda v: "%.3f" % v))
print("max |step diff| old %.3f unified %.3f" % (d.old_step_diff.abs().max(), d.unified_step_diff.abs().max()))
