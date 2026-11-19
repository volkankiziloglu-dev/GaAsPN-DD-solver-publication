"""Architecture-specific leakage, rebuilt with the Klitzke models of this project (dd_solver v1.1.2; 2026-11-13).
Replaces the inputs of analyze_leakage_archspecific_lux.py, analyze_leakage_r6.py and run_paper_leakage_archspecific.m, which used
the Klitzke_Validasyon models (RHJ N_D = 5e17, other optics):
  RHJ excess current dJ(V) = Jrec,meas - Jrec,model of the lower-doped Klitzke RHJ, model N_D = 2e16, mesh 4
      (results/paper_klitzke_rhj_excess_v112/rhj_excess_current_mesh4.csv; cases/run_klitzke_rhj_excess_v112.m), measured at
      100, 200, 300, 500 and 1000 lx (used directly; no interpolation in illuminance). Table points 0.80-1.20 V as before; dJ <= 0
      (model recombination at least as large as measured) is treated as no leakage; log-linear extrapolation outside 0.80-1.20 V
      from the two end points when both are positive, otherwise constant.
  HJ branch for HJ-0, F-HJ, RBSF-GaP: J0p*(exp(Vd/(n VT))-1), Vd = V - J*Rp, either the earlier parameters ('old': 4.97e-10 A/cm2,
      3.72, 8.82e3 Ohm cm2) or the refit on this project's HJ model ('new': results/paper_klitzke_rhj_excess_v112/hj_branch_fit.csv;
      cases/fit_klitzke_hj_branch_v112.m).
Curves: results/paper_voc_lux (600 and 10 ns; 100, 200, 500, 1000 lx) and results/paper_indoor_production (6500 K, 100 and 1000 lx).
Output: results/paper_leakage_v112/*.csv (old files untouched)."""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from unified_metrics import pv_metrics

root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EX = os.path.join(root, "results", "paper_klitzke_rhj_excess_v112")
D = os.path.join(root, "results", "paper_voc_lux"); P = os.path.join(root, "results", "paper_indoor_production")
OUT = os.path.join(root, "results", "paper_leakage_v112"); os.makedirs(OUT, exist_ok=True)
VT = 0.025852
X = pd.read_csv(os.path.join(EX, "rhj_excess_current_mesh4.csv"))
X = X[(X.V_V > 0.79) & (X.V_V < 1.21)].reset_index(drop=True)
Vt = X.V_V.values
DJ = {int(c.split("_")[1]): X[c].values * 1e-6 for c in X.columns if c.startswith("dJ_")}
F = pd.read_csv(os.path.join(EX, "hj_branch_fit.csv")).iloc[0]
BR = {"old": (4.97e-10, 3.72, 8.82e3), "new": (float(F.J0p_Acm2), float(F.n), float(F.Rp_Ohmcm2))}


def rhj_excess(V, lux):
    t = np.clip(DJ[lux], 1e-15, None); lnT = np.log(t)
    J = np.exp(np.interp(V, Vt, lnT))
    for (a, b, side) in ((0, 1, V < Vt[0]), (-2, -1, V > Vt[-1])):
        if DJ[lux][a] > 0 and DJ[lux][b] > 0:
            p = np.polyfit(Vt[[a, b]], lnT[[a, b]], 1); J = np.where(side, np.exp(np.polyval(p, V)), J)
    J[J <= 2e-15] = 0; J[V <= 0] = 0
    return J


def hj_branch(V, J0p, n, Rp):
    J = np.zeros_like(V, dtype=float)
    for i, v in enumerate(V):
        if v <= 0: continue
        lo, hi = 0.0, v
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            if mid + Rp * J0p * (np.exp(min(mid / (n * VT), 700)) - 1) - v > 0: hi = mid
            else: lo = mid
        J[i] = J0p * (np.exp(min(lo / (n * VT), 700)) - 1)
    return J


def leak(code, V, lux, sc, br):
    if sc == 0: return 0 * V
    if code == "RHJ": return sc * rhj_excess(V, lux)
    J0, n, Rp = BR[br]; return hj_branch(V, J0 * sc, n, Rp)


codes = ["HJ-0", "F-HJ", "RBSF-GaP", "RHJ"]
rows = []
for tau in (600, 10):
    for lux in (100, 200, 500, 1000):
        for code in codes:
            meta = pd.read_csv(os.path.join(D, f"voc_lux_{code}.csv")); q = meta[(meta.tauP_ns == tau) & (meta.lux == lux)]
            c = pd.read_csv(os.path.join(D, f"{code}_tau{tau}_{lux}lx.csv")); c = c[np.isfinite(c.V_V) & np.isfinite(c.J_Acm2)].sort_values("V_V")
            V, J, Pin = c.V_V.values, c.J_Acm2.values, float(q.Pin_Wm2.iloc[0])
            for br in ("old", "new"):
                for sc in (0, 0.5, 1, 2):
                    m = pv_metrics(V, J - leak(code, V, lux, sc, br), Pin)
                    rows.append(("voc_lux", tau, lux, code, br, sc, m["Jsc"] * 1e6, m["Voc"] * 1e3, m["FF"] * 100, m["PCE"]))
S = pd.read_csv(os.path.join(P, "indoor_production_summary.csv"))
for lux in (1000, 100):
    for code in codes:
        q = S[(S.Architecture == code) & (S.Spectrum == "6500K") & (S.Lux == lux)]; Pin = float(q.Pin_Wm2.iloc[0])
        c = pd.read_csv(os.path.join(P, f"{code}_6500K_{lux}lux_jv.csv")); c = c[np.isfinite(c.V_V) & np.isfinite(c.J_Acm2)].sort_values("V_V")
        V, J = c.V_V.values, c.J_Acm2.values
        for br in ("old", "new"):
            for sc in (0, 0.5, 1, 2):
                m = pv_metrics(V, J - leak(code, V, lux, sc, br), Pin)
                rows.append(("production", 600, lux, code, br, sc, m["Jsc"] * 1e6, m["Voc"] * 1e3, m["FF"] * 100, m["PCE"]))
T = pd.DataFrame(rows, columns=["Set", "tauP_ns", "Lux", "Architecture", "HJbranch", "Scale", "Jsc_uAcm2", "Voc_mV", "FF_pct", "PCE_pct"])
T.to_csv(os.path.join(OUT, "leakage_archspecific_v112.csv"), index=False, lineterminator=chr(10), float_format="%.5f")
W = T.pivot_table(index=["Set", "HJbranch", "tauP_ns", "Lux", "Scale"], columns="Architecture", values="PCE_pct").reset_index()
W["RBSF_minus_RHJ"] = W["RBSF-GaP"] - W["RHJ"]; W["best"] = W[codes].idxmax(axis=1)
W.to_csv(os.path.join(OUT, "leakage_archspecific_v112_wide.csv"), index=False, lineterminator=chr(10), float_format="%.3f")
pd.set_option("display.width", 220); pd.set_option("display.max_rows", 300)
print(W.round(2).to_string(index=False))
