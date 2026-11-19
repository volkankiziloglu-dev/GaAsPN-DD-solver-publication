"""Thickness- and mobility-normalised form of the design map.
Inputs: r*(tau_p) at absorber thickness W = 0.75, 1.5 and 2.25 um (mu_n = 500; results/paper_r4_analyses/DIM_*of8.csv and the nominal
merged design map) and at mu_n = 50, 100, 150, 500 (W = 1.5 um; results/paper_design_map_mu/rstar_mu.csv). mu_p = 112 cm2/Vs.
Dimensionless groups: electron diffusion length of the p-type absorber over W, L_e/W = sqrt(mu_n kT/q tau_p)/W; hole diffusion length of
the n-type absorber at the crossover, L_h*/W = sqrt(mu_p kT/q r* tau_p)/W; lifetime ratio at the crossover r* = tau_h*/tau_e.
r* is found as in analyze_design_map.py (linear interpolation of PCE_RHJ - PCE_RBSF in log10 r).
Output: results/paper_r4_analyses/dimensionless_rstar.csv and a printed summary."""
import glob, os
import numpy as np, pandas as pd

root = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); R = os.path.join(root, "results")
VT, MUP, W0 = 0.025852, 112.0, 1.5e-4


def rstar(rb, rh):
    rh = rh.sort_values("ratio"); r = rh.ratio.values; d = rh.PCE_pct.values - rb; ok = np.isfinite(d); r, d = r[ok], d[ok]
    if len(r) < 2: return np.nan
    if d[0] >= 0: return r[0] if d[0] == 0 else -r[0]          # negative: crossover below the smallest ratio of the grid
    if d[-1] < 0: return np.inf
    k = np.where(d >= 0)[0][0]; lr = np.log10(r)
    return 10 ** (lr[k - 1] + (0 - d[k - 1]) * (lr[k] - lr[k - 1]) / (d[k] - d[k - 1]))


rows = []
D = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(os.path.join(R, "paper_r4_analyses", "DIM_*of8.csv")))], ignore_index=True)
for (ts, cond, tp), g in D.groupby(["tScale", "Condition", "tauP_ns"]):
    rb = g[g.Architecture == "RBSF-GaP"]
    if rb.empty: continue
    rows.append(("W", 500.0, ts * 1.5, cond, tp, rstar(float(rb.PCE_pct.iloc[0]), g[g.Architecture == "RHJ"])))
M = pd.read_csv(os.path.join(R, "paper_design_map_mu", "rstar_mu.csv"))
for _, q in M.iterrows():
    rows.append(("W" if q.muN == 500 else "mu", float(q.muN), 1.5, q.Condition, q.tauP_ns, q.r_star))
T = pd.DataFrame(rows, columns=["series", "muN", "W_um", "Condition", "tauP_ns", "r_star"])
T["Le_W"] = np.sqrt(T.muN * VT * T.tauP_ns * 1e-9) / (T.W_um * 1e-4)
T["Lh_W"] = np.sqrt(MUP * VT * T.r_star.abs() * T.tauP_ns * 1e-9) / (T.W_um * 1e-4)
T["tauW_ns"] = (T.W_um * 1e-4) ** 2 / (MUP * VT) * 1e9          # hole transit (diffusion) time across W
T.to_csv(os.path.join(R, "paper_r4_analyses", "dimensionless_rstar.csv"), index=False, lineterminator=chr(10), float_format="%.4g")
pd.set_option("display.width", 200); pd.set_option("display.max_rows", 200)
for cond, g in T.groupby("Condition"):
    print("==", cond)
    print(g.pivot_table(index="tauP_ns", columns=["muN", "W_um"], values="r_star").round(2).to_string())
    print(g.sort_values("Le_W")[["muN", "W_um", "tauP_ns", "Le_W", "Lh_W", "r_star"]].round(3).to_string(index=False))
