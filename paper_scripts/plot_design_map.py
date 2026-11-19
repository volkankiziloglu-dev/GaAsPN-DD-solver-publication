"""Figure (with +-0.1-point band and mu_n = 100 curves): design map of the RHJ versus RBSF-GaP (crossover ratio r*(tau_p)) and the rear-layer gain over HJ-0.  Reads results/paper_design_map_merged/rstar.csv."""
import os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = pd.read_csv(os.path.join(root, "results", "paper_design_map_merged", "rstar.csv"))
BAND = pd.read_csv(os.path.join(root, "results", "paper_design_map_merged", "rstar_band.csv"))   # +-0.1-point indifference band (analyze_rstar_band.py)
MU = pd.read_csv(os.path.join(root, "results", "paper_design_map_mu", "rstar_mu.csv"))
plt.rcParams.update({"font.size": 9, "axes.linewidth": 0.8, "font.family": "DejaVu Sans"})
cA, cI = "#0072B2", "#D55E00"
fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.1), constrained_layout=True)
a = ax[0]
for cond, c, lab, mk in [("AM1.5G", cA, "AM1.5G", "o"), ("6500K_1000lx", cI, "6500 K, 1000 lx", "s")]:
    g = R[(R.Condition == cond) & (R.Sback == 1e5)].sort_values("tauP_ns")
    a.plot(g.tauP_ns, g.r_star, marker=mk, color=c, label=lab, lw=1.6, ms=4.5)
    B = BAND[BAND.Condition == cond].sort_values("tauP_ns"); lo = B.r_m01.where(B.r_m01 > 0, 1.0)
    a.fill_between(B.tauP_ns, lo, B.r_p01, color=c, alpha=0.18, lw=0)
    m = MU[(MU.Condition == cond) & (MU.muN == 100)].sort_values("tauP_ns")
    a.plot(m.tauP_ns, m.r_star, ls="--", color=c, lw=1.2, label=lab + r", $\mu_n$=100")
a.axhspan(1, 30, color="0.9", zorder=0); a.text(1.15, 22, "ratios sampled (1-30)", fontsize=7, color="0.35", va="center")
a.axhline(100, color="0.4", ls=":", lw=1); a.text(1.15, 120, "GaInP, low injection: about 100", fontsize=7, color="0.35", va="bottom")
a.set_xscale("log"); a.set_yscale("log"); a.set_ylim(0.8, 300); a.set_xlim(0.9, 700)
a.set_xlabel(r"p-type absorber lifetime $\tau_p$ (ns)"); a.set_ylabel(r"crossover ratio $r^*=\tau_n/\tau_p$")
a.text(200, 3.2, "RHJ ahead", fontsize=8, color="0.2", ha="center")
a.text(5, 1.12, "RBSF-GaP ahead", fontsize=8, color="0.2", ha="center", va="bottom")
a.legend(frameon=False, loc="upper right", fontsize=6.5, bbox_to_anchor=(1.0, 0.84)); a.set_title("(a) RHJ overtakes RBSF-GaP above $r^*$", fontsize=9, loc="left")
b = ax[1]
for cond, c, lab in [("AM1.5G", cA, "AM1.5G"), ("6500K_1000lx", cI, "1000 lx")]:
    for sb, ls, mk in [(1e5, "-", "o"), (1e3, "--", "^")]:
        g = R[(R.Condition == cond) & (R.Sback == sb)].sort_values("tauP_ns")
        b.plot(g.tauP_ns, g.PCE_RBSF - g.PCE_HJ0, ls=ls, marker=mk, color=c, lw=1.5, ms=4, label=f"{lab}, $S_{{back}}$=$10^{int(np.log10(sb))}$")
b.axhline(0, color="0.3", lw=0.8)
b.set_xscale("log"); b.set_xlabel(r"absorber lifetime $\tau_p$ (ns)"); b.set_ylabel("PCE(RBSF-GaP) $-$ PCE(HJ-0) (points)")
b.legend(frameon=False, fontsize=7, loc="upper left"); b.set_title("(b) rear-layer gain over the homojunction", fontsize=9, loc="left")
for x in ax: x.spines[["top", "right"]].set_visible(False)
out = os.path.join(root, "figures", "main_results"); os.makedirs(out, exist_ok=True)
fig.savefig(os.path.join(out, "fig_design_map.pdf")); fig.savefig(os.path.join(out, "fig_design_map.png"), dpi=200)
print("saved")
