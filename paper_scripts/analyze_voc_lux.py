"""Light-intensity dependence of V_oc (results/paper_voc_lux): slope dV_oc/dlog10(lux) and ideality n = slope / (kT/q ln 10) for the four architectures at two lifetimes.
Extractor: cases/unified_metrics.py on all finite points. Slopes: two-point 100 -> 1000 lx (the interval of the Klitzke measurement) and a least-squares fit over 100-5000 lx.
Output: results/paper_voc_lux/voc_lux_summary.csv, voc_lux_slopes.csv"""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from unified_metrics import pv_metrics
root = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); d = os.path.join(root, "results", "paper_voc_lux")
kTq = 8.617333262e-5 * 300; per_dec = kTq * np.log(10) * 1e3   # mV per decade for n = 1
rows = []
for code in ["HJ-0", "F-HJ", "RBSF-GaP", "RHJ"]:
    meta = pd.read_csv(os.path.join(d, f"voc_lux_{code}.csv"))
    for _, r in meta.iterrows():
        f = os.path.join(d, f"{code}_tau{r.tauP_ns:g}_{r.lux:g}lx.csv")
        if not os.path.exists(f): continue
        c = pd.read_csv(f); c = c[np.isfinite(c.J_Acm2)]
        m = pv_metrics(c.V_V, c.J_Acm2, r.Pin_Wm2)
        rows.append((code, r.tauP_ns, r.lux, m["Jsc"] * 1e6, m["Voc"] * 1e3, m["FF"] * 100, m["PCE"], r.Class))
S = pd.DataFrame(rows, columns=["Architecture", "tauP_ns", "lux", "Jsc_uAcm2", "Voc_mV", "FF_pct", "PCE_pct", "Class"])
S.to_csv(os.path.join(d, "voc_lux_summary.csv"), index=False, float_format="%.4f")
sl = []
for (a, t), g in S.groupby(["Architecture", "tauP_ns"]):
    g = g.sort_values("lux"); x = np.log10(g.lux.values); y = g.Voc_mV.values
    def two(l1, l2):
        y1 = g[g.lux == l1].Voc_mV.iloc[0]; y2 = g[g.lux == l2].Voc_mV.iloc[0]; return (y2 - y1) / np.log10(l2 / l1)
    ok = (g.lux >= 100).values
    fit = np.polyfit(x[ok], y[ok], 1)[0]
    # dJsc exponent (Jsc ~ lux^alpha)
    alpha = np.polyfit(np.log10(g.lux.values[ok]), np.log10(g.Jsc_uAcm2.values[ok]), 1)[0]
    sl.append((a, t, two(100, 1000), two(100, 1000) / per_dec, fit, fit / per_dec, alpha, two(1000, 5000), g.Class.iloc[(g.lux == 1000).values.argmax()]))
L = pd.DataFrame(sl, columns=["Architecture", "tauP_ns", "slope_100_1000", "n_100_1000", "slope_fit_100_5000", "n_fit", "Jsc_exponent", "slope_1000_5000", "Class_1000lx"])
L.to_csv(os.path.join(d, "voc_lux_slopes.csv"), index=False, float_format="%.3f")
pd.set_option("display.width", 220)
print("n = 1 corresponds to %.1f mV/dec" % per_dec)
print(S.pivot_table(index=["Architecture", "tauP_ns"], columns="lux", values="Voc_mV").round(1).to_string())
print(L.to_string(index=False, float_format=lambda v: "%.2f" % v))
