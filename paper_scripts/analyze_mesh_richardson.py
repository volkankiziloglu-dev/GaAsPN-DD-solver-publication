"""Three-level mesh study with one common extractor (cases/unified_metrics.py) and Richardson extrapolation.
Reads results/paper_mesh_richardson/<code>_x<scale>_<cond>.csv (all finite points), computes J_sc, V_oc, FF, PCE for scale 1, 2, 4 (nodes), the observed order
p = log2((f1 - f2)/(f2 - f4)) and the extrapolated value f_inf = f4 + (f4 - f2)/(2^p - 1) when 0.5 < p < 3 (otherwise the sequence is not asymptotic and only the differences are given).
Outputs: results/paper_mesh_richardson/richardson_summary.csv and richardson_table.tex (rows for the supplement)."""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from unified_metrics import pv_metrics

root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
d = os.path.join(root, "results", "paper_mesh_richardson")
Pin = {"AM1.5G": None, "6500K_1000lx": None}
sys.path.insert(0, os.path.join(root, "paper_scripts"))
# input powers from the production summaries
S = pd.read_csv(os.path.join(root, "results", "paper_indoor_production", "indoor_production_summary.csv"))
Pin["6500K_1000lx"] = float(S[(S.Architecture == "RBSF-GaP") & (S.Spectrum == "6500K") & (S.Lux == 1000)].Pin_Wm2.iloc[0])
Pin["AM1.5G"] = 1000.3706555734421
rows = []
for code in ["HJ-0", "F-HJ", "RBSF-GaP", "RHJ"]:
    for cond in ["AM1.5G", "6500K_1000lx"]:
        mm = {}
        for sc in (1, 2, 4):
            f = os.path.join(d, f"{code}_x{sc}_{cond}.csv")
            if not os.path.exists(f): continue
            c = pd.read_csv(f); c = c[np.isfinite(c.J_Acm2)]
            mm[sc] = pv_metrics(c.V_V, c.J_Acm2, Pin[cond])
        if len(mm) < 3: continue
        for q, scale in (("Jsc", 1e3), ("Voc", 1e3), ("FF", 100), ("PCE", 1)):
            f1, f2, f4 = (mm[s][q] * scale for s in (1, 2, 4))
            den = f2 - f4
            p = np.log2((f1 - f2) / den) if den != 0 and (f1 - f2) / den > 0 else np.nan
            finf = f4 + (f4 - f2) / (2 ** p - 1) if np.isfinite(p) and 0.5 < p < 3 else np.nan
            rows.append((code, cond, q, f1, f2, f4, f1 - f2, f2 - f4, p, finf, (f1 - finf) if np.isfinite(finf) else np.nan))
R = pd.DataFrame(rows, columns=["Architecture", "Condition", "metric", "nominal", "x2", "x4", "d12", "d24", "order", "extrapolated", "nominal_minus_extrapolated"])
R.to_csv(os.path.join(d, "richardson_summary.csv"), index=False, float_format="%.5f")
pd.set_option("display.width", 220)
print(R.to_string(index=False, float_format=lambda v: "%.4f" % v))
