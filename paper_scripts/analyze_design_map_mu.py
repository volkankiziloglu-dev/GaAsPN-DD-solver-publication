"""Sensitivity of the crossover ratio r*(tau_p) to the assumed absorber electron mobility mu_n (50, 100, 150 and the nominal 500 cm2/Vs; mu_p = 112 cm2/Vs).
Reads results/paper_design_map_mu/design_map_mu<mu>_{a,b}.csv (RBSF-GaP and RHJ, S_back = 1e5; AM1.5G and 6500 K/1000 lx) and, for mu_n = 500, the merged design map.
r* is found as in analyze_design_map.py (linear interpolation of PCE_RHJ - PCE_RBSF versus log10 r; r* < 1: RHJ ahead at equal lifetimes; > 60: no crossover in the grid).
Output: results/paper_design_map_mu/rstar_mu.csv and a printed table."""
import glob, os
import numpy as np, pandas as pd
root = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); R = os.path.join(root, "results")
frames = []
for mu in (50, 100, 150):
    fs = sorted(glob.glob(os.path.join(R, "paper_design_map_mu", f"design_map_mu{mu}_[ab].csv")))
    if fs:
        d = pd.concat([pd.read_csv(f) for f in fs]); d["muN"] = mu; frames.append(d)
m500 = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(os.path.join(R, "paper_design_map_merged", "design_map_[abcd].csv")))])
m500 = m500[(m500.Sback == 1e5) & m500.Architecture.isin(["RBSF-GaP", "RHJ"])].copy(); m500["muN"] = 500; frames.append(m500)
D = pd.concat(frames, ignore_index=True)
rows = []
for (mu, cond, tp), g in D.groupby(["muN", "Condition", "tauP_ns"]):
    rb = g[g.Architecture == "RBSF-GaP"]
    if rb.empty: continue
    rb = rb.iloc[0]; rh = g[g.Architecture == "RHJ"].sort_values("ratio")
    r = rh.ratio.values; diff = (rh.PCE_pct - rb.PCE_pct).values; ok = np.isfinite(diff); r, diff = r[ok], diff[ok]
    if len(r) < 2: continue
    lr = np.log10(r)
    if diff[0] >= 0: rstar = 1.0 if diff[0] == 0 else np.nan
    elif diff[-1] < 0: rstar = np.inf
    else:
        k = np.where(diff >= 0)[0][0]; rstar = 10 ** (lr[k - 1] + (0 - diff[k - 1]) * (lr[k] - lr[k - 1]) / (diff[k] - diff[k - 1]))
    rows.append((mu, cond, tp, rb.PCE_pct, diff[0], rstar))
R2 = pd.DataFrame(rows, columns=["muN", "Condition", "tauP_ns", "PCE_RBSF", "RHJ_minus_RBSF_r1", "r_star"])
R2.to_csv(os.path.join(R, "paper_design_map_mu", "rstar_mu.csv"), index=False, float_format="%.3f")
pd.set_option("display.width", 200)
for cond in ("AM1.5G", "6500K_1000lx"):
    p = R2[R2.Condition == cond].pivot(index="tauP_ns", columns="muN", values="r_star")
    print(cond, "r*"); print(p.round(2).to_string())
    q = R2[R2.Condition == cond].pivot(index="tauP_ns", columns="muN", values="RHJ_minus_RBSF_r1")
    print(cond, "RHJ - RBSF at r = 1 (points)"); print(q.round(2).to_string())
