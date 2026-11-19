"""Perimeter recombination and loss-transferred efficiency estimate.

1. Perimeter diode. A mesa perimeter adds J_per(V) = J0_per * (P/A) * (exp(V/(2 V_T)) - 1) (A/cm2; J0_per in A/cm of perimeter,
   ideality 2) to the dark current. It is superposed on the model J-V curves (nominal parameters, tau = 0.6 us, CIE 6500 K;
   results/paper_voc_lux/<code>_tau600_<lux>lx.csv) for 50-5000 lx and an effective J02 = J0_per * P/A swept over 1e-16 - 1e-9 A/cm2.
   For a square cell P/A = 4/sqrt(A): 40 cm-1 (1 mm2), 4 cm-1 (1 cm2). Output: PCE loss per architecture and the effective J02
   (and J0_per for 1 mm2 / 1 cm2) at which the PCE drops by 0.5 and 1 point; RBSF-GaP minus the others with the same J02
   (identical perimeter for all architectures).
2. Loss-transferred estimate. The measured/model PCE ratio of the Klitzke GaInP cells at 1000 and 100 lx (model without terminal
   losses, Table tab:klitzke-validation: HJ 34.5/38.90 and 24.7/35.27, RHJ 41.0/47.64 and 37.5/44.22) applied to the GaAsPN
   model PCE (homojunction ratio to HJ-0, F-HJ, RBSF-GaP; RHJ ratio to RHJ; both ratios given for RBSF-GaP); 100 lx uses the
   nearest model illuminance of the series (100 lx).
Output: results/paper_perimeter/perimeter_loss.csv, perimeter_thresholds.csv, loss_transferred.csv (printed)."""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from unified_metrics import pv_metrics

root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D = os.path.join(root, "results", "paper_voc_lux"); OUT = os.path.join(root, "results", "paper_perimeter"); os.makedirs(OUT, exist_ok=True)
VT = 0.025852
A = ["HJ-0", "F-HJ", "RBSF-GaP", "RHJ"]
LUX = [50, 100, 200, 500, 1000, 2000, 5000]
P = pd.read_csv(os.path.join(D, "lux_series_pce.csv")).drop_duplicates(["Architecture", "lux"])
Pin = {int(r.lux): r.Pin_Wm2 for r in P[P.Architecture == "RBSF-GaP"].itertuples()}
J02s = np.logspace(-16, -9, 141)
rows = []
for a in A:
    for lx in LUX:
        c = pd.read_csv(os.path.join(D, f"{a}_tau600_{lx}lx.csv")); c = c[np.isfinite(c.J_Acm2)]
        V = c.V_V.values; J = c.J_Acm2.values
        m0 = pv_metrics(V, J, Pin[lx])
        for j02 in J02s:
            m = pv_metrics(V, J - j02 * (np.exp(V / (2 * VT)) - 1), Pin[lx])
            rows.append((a, lx, j02, m0["PCE"], m["PCE"], m["Voc"] * 1e3, m["FF"] * 100))
L = pd.DataFrame(rows, columns=["Architecture", "lux", "J02_eff_Acm2", "PCE0", "PCE", "Voc_mV", "FF_pct"])
L["loss"] = L.PCE0 - L.PCE
L.to_csv(os.path.join(OUT, "perimeter_loss.csv"), index=False, float_format="%.6g")


def thr(g, d):
    x = np.log10(g.J02_eff_Acm2.values); y = g.loss.values
    k = np.where(y >= d)[0]
    if len(k) == 0: return np.nan
    k = k[0]
    if k == 0: return 10 ** x[0]
    return 10 ** (x[k - 1] + (d - y[k - 1]) * (x[k] - x[k - 1]) / (y[k] - y[k - 1]))


T = []
for (a, lx), g in L.groupby(["Architecture", "lux"]):
    for d in (0.5, 1.0):
        j = thr(g, d); T.append((a, lx, d, j, j / 40, j / 4))
T = pd.DataFrame(T, columns=["Architecture", "lux", "PCE_loss_pts", "J02_eff_Acm2", "J0per_1mm2_Acm", "J0per_1cm2_Acm"])
T.to_csv(os.path.join(OUT, "perimeter_thresholds.csv"), index=False, float_format="%.3g")
pd.set_option("display.width", 220)
print("effective J02 (A/cm2) for a 1-point PCE loss:"); print(T[T.PCE_loss_pts == 1.0].pivot(index="lux", columns="Architecture", values="J02_eff_Acm2").to_string(float_format="%.2e"))
print("effective J02 (A/cm2) for a 0.5-point PCE loss:"); print(T[T.PCE_loss_pts == 0.5].pivot(index="lux", columns="Architecture", values="J02_eff_Acm2").to_string(float_format="%.2e"))
# differences between architectures at a common J02
for j in (1e-16, 3e-16, 1e-15, 3e-15):
    k = np.argmin(np.abs(np.log10(J02s) - np.log10(j))); jj = J02s[k]
    s = L[np.isclose(L.J02_eff_Acm2, jj, rtol=1e-9, atol=0)].pivot(index="lux", columns="Architecture", values="PCE")
    print("J02_eff = %.0e A/cm2: PCE" % jj); print(s.round(2).to_string())
# 2. loss-transferred estimate
ratio = {("HJ", 1000): 34.5 / 38.90, ("HJ", 100): 24.7 / 35.27, ("RHJ", 1000): 41.0 / 47.64, ("RHJ", 100): 37.5 / 44.22}
E = []
for lx in (1000, 100):
    for a in A:
        p0 = float(P[(P.Architecture == a) & (P.lux == lx)].PCE_pct.iloc[0])
        for src in ("HJ", "RHJ"):
            E.append((lx, a, p0, src, ratio[(src, lx)], p0 * ratio[(src, lx)]))
E = pd.DataFrame(E, columns=["lux", "Architecture", "PCE_model", "ratio_source", "ratio", "PCE_loss_transferred"])
E.to_csv(os.path.join(OUT, "loss_transferred.csv"), index=False, float_format="%.3f")
print(E.round(3).to_string())
