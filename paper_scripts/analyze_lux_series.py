"""PCE versus illuminance (50-5000 lx, CIE 6500 K) of the four architectures from results/paper_voc_lux, with and without an unmodelled leakage path.
Leakage cases (as in run_paper_leakage_robustness.m, post-processing of the intrinsic curves): ohmic shunt Rsh = 1e7, 1e6, 1e5 Ohm cm2, and the diode-plus-series-resistor branch fitted to the excess current of the
measured Klitzke homojunction (J0p = 4.97e-10 A/cm2, np = 3.72, Rp = 8.82e3 Ohm cm2). Curves are pchip-interpolated to 2 mV before the subtraction; metrics with the common extractor.
Output: results/paper_voc_lux/lux_series_pce.csv, lux_series_leakage.csv"""
import os, sys
import numpy as np, pandas as pd
from scipy.interpolate import PchipInterpolator
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from unified_metrics import pv_metrics
root = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); d = os.path.join(root, "results", "paper_voc_lux")
VT = 0.025852


def leak_branch(V, J0p, np_, Rp):
    J = np.zeros_like(V)
    for i, v in enumerate(V):
        if v <= 0: continue
        lo, hi = 0.0, v
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            f = mid + Rp * J0p * (np.exp(min(mid / (np_ * VT), 700)) - 1) - v
            if f > 0: hi = mid
            else: lo = mid
        J[i] = J0p * (np.exp(min(lo / (np_ * VT), 700)) - 1)
    return J


cases = [("intrinsic", None), ("Rsh 1e7", 1e7), ("Rsh 1e6", 1e6), ("Rsh 1e5", 1e5), ("branch", "b")]
rows = []; lrows = []
for code in ["HJ-0", "F-HJ", "RBSF-GaP", "RHJ"]:
    meta = pd.read_csv(os.path.join(d, f"voc_lux_{code}.csv"))
    for _, r in meta[meta.tauP_ns == 600].iterrows():
        f = os.path.join(d, f"{code}_tau600_{r.lux:g}lx.csv")
        if not os.path.exists(f): continue
        c = pd.read_csv(f); c = c[np.isfinite(c.J_Acm2)].sort_values("V_V")
        m = pv_metrics(c.V_V, c.J_Acm2, r.Pin_Wm2)
        rows.append((code, r.lux, r.Pin_Wm2, m["Jsc"] * 1e6, m["Voc"] * 1e3, m["FF"] * 100, m["PCE"], r.Class))
        V = np.arange(0, c.V_V.max() + 1e-9, 0.002); J = PchipInterpolator(c.V_V.values, c.J_Acm2.values)(V)
        for nm, p in cases:
            Jl = np.zeros_like(V)
            if isinstance(p, float): Jl = V / p
            elif p == "b": Jl = leak_branch(V, 4.97e-10, 3.72, 8.82e3)
            mm = pv_metrics(V, J - Jl, r.Pin_Wm2)
            lrows.append((code, r.lux, nm, mm["Jsc"] * 1e6, mm["Voc"] * 1e3, mm["FF"] * 100, mm["PCE"], r.Class))
S = pd.DataFrame(rows, columns=["Architecture", "lux", "Pin_Wm2", "Jsc_uAcm2", "Voc_mV", "FF_pct", "PCE_pct", "Class"])
S.to_csv(os.path.join(d, "lux_series_pce.csv"), index=False, float_format="%.4f")
L = pd.DataFrame(lrows, columns=["Architecture", "lux", "case", "Jsc_uAcm2", "Voc_mV", "FF_pct", "PCE_pct", "Class"])
L.to_csv(os.path.join(d, "lux_series_leakage.csv"), index=False, float_format="%.4f")
pd.set_option("display.width", 220)
print(S.pivot_table(index="lux", columns="Architecture", values="PCE_pct").round(2).to_string())
print(S.pivot_table(index="lux", columns="Architecture", values="Class", aggfunc="first").to_string())
for lux in (200, 500, 1000):
    p = L[L.lux == lux].pivot_table(index="case", columns="Architecture", values="PCE_pct").loc[[c[0] for c in cases]]
    p["RBSF-RHJ"] = p["RBSF-GaP"] - p["RHJ"]; p["RHJ-F-HJ"] = p["RHJ"] - p["F-HJ"]; p["RHJ-HJ0"] = p["RHJ"] - p["HJ-0"]
    print("\n== %d lx ==" % lux); print(p.round(2).to_string())
