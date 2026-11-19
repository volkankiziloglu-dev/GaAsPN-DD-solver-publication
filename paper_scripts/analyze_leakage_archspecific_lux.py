"""Architecture-specific leakage at 100, 200, 500 and 1000 lx (CIE 6500 K), round-5 review R3-1.
Same post-processing as cases/run_paper_leakage_archspecific.m, applied to the lux-series curves (results/paper_voc_lux, nominal 600 ns and
10 ns, native 0.02 V grid, unified extractor):
  HJ-type designs (HJ-0, F-HJ, RBSF-GaP): diode + series-resistor branch fitted to the excess current of the Klitzke homojunction
  (J0p = 4.97e-10 A/cm2 x scale, n = 3.72, Rp = 8.82e3 Ohm cm2; light-independent);
  RHJ: the measured excess current dJ(V) of the lower-doped Klitzke RHJ, tabulated at 100 and 1000 lx; at 200 and 500 lx it is interpolated
  per voltage as a power law in illuminance, dJ(L) = dJ100 * (L/100)^k(V), k = log10(dJ1000/dJ100) (an assumption: no measurement exists).
Modes: archspecific, rhj_branch_all, hj_branch_all; scales 0 (intrinsic), 0.5, 1, 2.
Output: results/paper_leakage_robustness/leakage_archspecific_lux.csv and a ranking summary."""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from unified_metrics import pv_metrics

root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D = os.path.join(root, "results", "paper_voc_lux"); OUT = os.path.join(root, "results", "paper_leakage_robustness")
VT = 0.025852
Vt = np.array([0.80, 0.90, 1.00, 1.10, 1.20])
dJ = {100: np.array([0.12, 0.23, 0.38, 0.61, 2.00]) * 1e-6, 1000: np.array([0.99, 1.32, 1.69, 2.93, 5.91]) * 1e-6}


def rhj_table(lux):
    if lux in dJ: return dJ[lux]
    k = np.log10(dJ[1000] / dJ[100])
    return dJ[100] * (lux / 100.0) ** k


def rhj_excess(V, lux):
    lnT = np.log(rhj_table(lux))
    p1 = np.polyfit(Vt[:2], lnT[:2], 1); p2 = np.polyfit(Vt[3:], lnT[3:], 1)
    lnJ = np.interp(V, Vt, lnT)
    lnJ = np.where(V < Vt[0], np.polyval(p1, V), lnJ); lnJ = np.where(V > Vt[-1], np.polyval(p2, V), lnJ)
    J = np.exp(lnJ); J[V <= 0] = 0
    return J


def leak_branch(V, J0p, np_, Rp):
    J = np.zeros_like(V)
    for i, v in enumerate(V):
        if v <= 0: continue
        lo, hi = 0.0, v
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            f = mid + Rp * J0p * (np.exp(min(mid / (np_ * VT), 700)) - 1) - v
            if f > 0: hi = mid
            else: lo = mid
        J[i] = J0p * (np.exp(min(lo / (np_ * VT), 700)) - 1)
    return J


codes = ["HJ-0", "F-HJ", "RBSF-GaP", "RHJ"]
rows = []
for tau in (600, 10):
    for code in codes:
        meta = pd.read_csv(os.path.join(D, f"voc_lux_{code}.csv"))
        for lux in (100, 200, 500, 1000):
            q = meta[(meta.tauP_ns == tau) & (meta.lux == lux)]
            f = os.path.join(D, f"{code}_tau{tau}_{lux}lx.csv")
            if q.empty or not os.path.exists(f): continue
            Pin = float(q.Pin_Wm2.iloc[0]); cls = q.Class.iloc[0]
            c = pd.read_csv(f); c = c[np.isfinite(c.V_V) & np.isfinite(c.J_Acm2)].sort_values("V_V")
            V = c.V_V.values; J = c.J_Acm2.values
            for sc in (0, 0.5, 1, 2):
                for mode in ("archspecific", "rhj_branch_all", "hj_branch_all"):
                    if sc == 0 and mode != "archspecific": continue
                    useRHJ = mode == "rhj_branch_all" or (mode == "archspecific" and code == "RHJ")
                    if sc == 0: Jl = np.zeros_like(V)
                    elif useRHJ: Jl = sc * rhj_excess(V, lux)
                    else: Jl = leak_branch(V, 4.97e-10 * sc, 3.72, 8.82e3)
                    m = pv_metrics(V, J - Jl, Pin)
                    rows.append((tau, code, lux, cls, mode, sc, m["Jsc"] * 1e6, m["Voc"] * 1e3, m["FF"] * 100, m["PCE"]))
T = pd.DataFrame(rows, columns=["tauP_ns", "Architecture", "Lux", "Class", "Mode", "Scale", "Jsc_uAcm2", "Voc_mV", "FF_pct", "PCE_pct"])
T.to_csv(os.path.join(OUT, "leakage_archspecific_lux.csv"), index=False, lineterminator=chr(10))
for tau in (600, 10):
    for lux in (100, 200, 500, 1000):
        for mode in ("archspecific", "rhj_branch_all", "hj_branch_all"):
            for sc in (0, 0.5, 1, 2):
                if sc == 0 and mode != "archspecific": continue
                g = T[(T.tauP_ns == tau) & (T.Lux == lux) & (T.Mode == mode) & (T.Scale == sc)].set_index("Architecture").PCE_pct
                if len(g) < 4: continue
                order = ">".join(g.sort_values(ascending=False).index)
                print("tau %3d %4d lx %-15s x%-3g HJ-0 %6.2f F-HJ %6.2f RBSF %6.2f RHJ %6.2f | %s | RBSF-RHJ %+.2f" % (
                    tau, lux, mode, sc, g["HJ-0"], g["F-HJ"], g["RBSF-GaP"], g["RHJ"], order, g["RBSF-GaP"] - g["RHJ"]))
