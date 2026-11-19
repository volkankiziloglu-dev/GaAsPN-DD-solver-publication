"""LaTeX tables of the group-B robustness analyses, read from results/paper_r4_analyses/ (cases/analyze_r4_groupB.py)
and results/paper_perimeter/ (cases/analyze_perimeter_loss.py). Writes latex/tables/supp_r4B_<name>.tex."""
import os
import numpy as np, pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
A = os.path.join(ROOT, "results", "paper_r4_analyses"); P = os.path.join(ROOT, "results", "paper_perimeter"); OUT = os.path.join(ROOT, "latex", "tables")
TEX = {"HJ-0": "HJ-0", "HJ-U": "HJ-U", "F-HJ": "F-HJ", "RBSF-GaP": "RBSF--GaP", "RHJ": "RHJ"}
CN = {"AM1.5G": "AM1.5G", "6500K_1000lx": "1000~lx", "6500K_200lx": "200~lx"}


def write(name, L):
    open(os.path.join(OUT, "supp_r4B_%s.tex" % name), "w", encoding="utf-8", newline="\n").write("\n".join(L) + "\n")


def f(x, nd=2):
    if not np.isfinite(x): return "--"
    s = ("%." + str(nd) + "f") % abs(x)
    return ("$-$" + s) if x < 0 and float(s) != 0 else s


# ---- 200 lx: r* and rear-contact controls ----
r = pd.read_csv(os.path.join(A, "rstar_200lx.csv")); sb = pd.read_csv(os.path.join(A, "sback_200lx.csv"))
L = [r"\begin{table}[htbp]", r"\centering",
     r"\caption{CIE 6500~K, 200~lx: crossover ratio $r^*(\tau_p)$ (RHJ versus RBSF--GaP, $S_{\rm back}=10^5$~cm\,s$^{-1}$) compared with 1000~lx and AM1.5G (design map), and the PCE difference (points) of RBSF--GaP to HJ-0, HJ-U and F-HJ at $S_{\rm back}=10^5$ / $10^3$~cm\,s$^{-1}$ (equal lifetimes).}",
     r"\label{tab:supp-r4-200lx}", r"\scriptsize", r"\begin{tabular}{rccccccc}", r"\toprule",
     r"$\tau_p$ (ns) & $r^*$ 200~lx & $r^*$ 1000~lx & $r^*$ AM1.5G & PCE RBSF--GaP (\%) & $-$HJ-0 & $-$HJ-U & $-$F-HJ \\", r"\midrule"]
for t in r.itertuples():
    a = sb[(sb.tauP_ns == t.tauP_ns) & (sb.Sback == 1e5)].iloc[0]; b = sb[(sb.tauP_ns == t.tauP_ns) & (sb.Sback == 1e3)].iloc[0]
    L.append("%g & %.2f & %.2f & %.2f & %.2f & %s / %s & %s / %s & %s / %s \\\\" % (t.tauP_ns, t.r_star_200lx, t.r_star_1000lx, t.r_star_AM, t.PCE_RBSF_200lx,
             f(a["RBSF-HJ-0"]), f(b["RBSF-HJ-0"]), f(a["RBSF-HJ-U"]), f(b["RBSF-HJ-U"]), f(a["RBSF-F-HJ"]), f(b["RBSF-F-HJ"])))
L += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]; write("200lx", L)

# ---- thickness / doping optimum ----
o = pd.read_csv(os.path.join(A, "opt_summary.csv"))
L = [r"\begin{table}[htbp]", r"\centering",
     r"\caption{Per-architecture optimum on a thickness $\times$ doping grid (absorber/base 0.375--3.0~$\mu$m, $10^{16}$--$10^{19}$~cm$^{-3}$; 38 points per architecture) at equal absorber lifetimes of 1 and 10~ns ($S_{\rm back}=10^5$~cm\,s$^{-1}$, RHJ $r=1$): PCE of the nominal stack (1.5~$\mu$m, $10^{17}$~cm$^{-3}$), best PCE, its thickness and doping. Lifetime and mobility are doping-independent and band-gap narrowing is not modelled, so the gains at $\geq10^{18}$~cm$^{-3}$ are upper bounds.}",
     r"\label{tab:supp-r4-opt}", r"\scriptsize", r"\begin{tabular}{llcccccc}", r"\toprule",
     r"$\tau$ & Architecture & \multicolumn{3}{c}{AM1.5G} & \multicolumn{3}{c}{1000~lx} \\",
     r" & & nominal & best & ($\mu$m, cm$^{-3}$) & nominal & best & ($\mu$m, cm$^{-3}$) \\", r"\midrule"]
for tp in (1, 10):
    for a in ["HJ-0", "HJ-U", "F-HJ", "RBSF-GaP", "RHJ"]:
        cells = []
        for c in ("AM1.5G", "6500K_1000lx"):
            x = o[(o.tauP_ns == tp) & (o.Condition == c) & (o.Architecture == a)].iloc[0]
            e = int(np.floor(np.log10(x.best_Nabs))); m = x.best_Nabs / 10 ** e
            cells += ["%.2f" % x.PCE_nominal_stack, "%.2f" % x.PCE_best, "%.2f, %s$10^{%d}$" % (x.best_thickness_um, "" if abs(m - 1) < 1e-9 else "%g$\\times$" % m, e)]
        L.append("%s & %s & %s \\\\" % ("%g~ns" % tp if a == "HJ-0" else "", TEX[a], " & ".join(cells)))
    if tp == 1: L.append(r"\midrule")
L += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]; write("opt", L)

# ---- single SRH centre ----
s2 = pd.read_csv(os.path.join(A, "srh2_rstar.csv")); s1 = pd.read_csv(os.path.join(A, "srh2_summary.csv"))
L = [r"\begin{table}[htbp]", r"\centering",
     r"\caption{Same SRH centre in the p-type (RBSF--GaP) and n-type (RHJ) absorber: crossover in $r=\tau_{p0}/\tau_{n0}$ ($E_t=E_i$) compared with the design map at $\tau_p=\tau_{n0}$, and RHJ$-$RBSF--GaP (points) at $r=10$ for $E_t-E_i=-0.8$, 0 and $+0.8$~eV.}",
     r"\label{tab:supp-r4-srh}", r"\scriptsize", r"\begin{tabular}{rlccccc}", r"\toprule",
     r"$\tau_{n0}$ (ns) & Condition & $r^*$ single centre & $r^*$ map & \multicolumn{3}{c}{RHJ$-$RBSF--GaP at $r=10$} \\",
     r" & & & & $-0.8$~eV & 0 & $+0.8$~eV \\", r"\midrule"]
for t in s2.itertuples():
    g = s1[(s1.tauP_ns == t.tau_n0_ns) & (s1.ratio == 10) & (s1.Condition == t.Condition)].set_index("EtRelEi_eV").RHJ_minus_RBSF
    L.append("%g & %s & %.1f & %.2f & %s & %s & %s \\\\" % (t.tau_n0_ns, CN[t.Condition], t.r_star_single_centre, t.r_star_map, f(g.get(-0.8, np.nan)), f(g.get(0.0, np.nan)), f(g.get(0.8, np.nan))))
L += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]; write("srh", L)

# ---- perimeter and loss-transferred ----
th = pd.read_csv(os.path.join(P, "perimeter_thresholds.csv")); lt = pd.read_csv(os.path.join(P, "loss_transferred.csv"))
L = [r"\begin{table}[htbp]", r"\centering",
     r"\caption{Effective perimeter saturation current $J_{02}=J_{0,\rm per}P/A$ (A\,cm$^{-2}$, ideality 2) that lowers the PCE by 1 point (nominal parameters, 6500~K); for a square cell $J_{0,\rm per}=J_{02}/40$ (1~mm$^2$) or $J_{02}/4$~A\,cm$^{-1}$ (1~cm$^2$).}",
     r"\label{tab:supp-r4-perimeter}", r"\scriptsize", r"\begin{tabular}{lcccc}", r"\toprule", r" & HJ-0 & F-HJ & RBSF--GaP & RHJ \\", r"\midrule"]
for lx in (50, 200, 1000, 5000):
    g = th[(th.lux == lx) & (th.PCE_loss_pts == 1.0)].set_index("Architecture").J02_eff_Acm2
    L.append("%d~lx & %s \\\\" % (lx, " & ".join(("%.1f$\\times10^{%d}$" % (g[a] / 10 ** np.floor(np.log10(g[a])), np.floor(np.log10(g[a])))) for a in ["HJ-0", "F-HJ", "RBSF-GaP", "RHJ"])))
# loss-transferred rows removed 2026-11-08 (estimate dropped from the paper; data remain in results/paper_perimeter/loss_transferred.csv)
L +=[r"\bottomrule", r"\end{tabular}", r"\end{table}"]; write("perimeter", L)
print("written")
