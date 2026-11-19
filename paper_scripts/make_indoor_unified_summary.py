"""Indoor production summary with the unified extractor (cases/unified_metrics.py) applied to the saved 0.02 V production curves (all finite points).
Reads results/paper_indoor_production/indoor_production_summary.csv and *_jv.csv; writes indoor_production_summary_unified.csv with the same columns (Jsc, Voc, FF, PCE, Pmax replaced)."""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from unified_metrics import pv_metrics
root = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); d = os.path.join(root, "results", "paper_indoor_production")
S = pd.read_csv(os.path.join(d, "indoor_production_summary.csv"))
U = S.copy()
for i, r in S.iterrows():
    f = os.path.join(d, f"{r.Architecture}_{r.Spectrum}_{int(r.Lux)}lux_jv.csv")
    if not os.path.exists(f): continue
    c = pd.read_csv(f); c = c[np.isfinite(c.J_Acm2)]
    m = pv_metrics(c.V_V, c.J_Acm2, r.Pin_Wm2)
    U.loc[i, ["Jsc_uAcm2", "Voc_mV", "FF_pct", "PCE_pct", "Pmax_mWcm2"]] = [m["Jsc"] * 1e6, m["Voc"] * 1e3, m["FF"] * 100, m["PCE"], m["Pmax"] * 1e3]
U.to_csv(os.path.join(d, "indoor_production_summary_unified.csv"), index=False, float_format="%.6g")
a = S[S.Architecture.isin(["HJ-0", "F-HJ", "RBSF-GaP", "RHJ"]) & (S.Lux == 1000)][["Architecture", "Spectrum", "PCE_pct"]].merge(
    U[["Architecture", "Spectrum", "Lux", "PCE_pct"]][U.Lux == 1000], on=["Architecture", "Spectrum"], suffixes=("_old", "_unified"))
print(a.to_string(index=False, float_format=lambda v: "%.3f" % v))
