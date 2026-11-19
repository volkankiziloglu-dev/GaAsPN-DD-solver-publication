"""Supplementary table: crossover ratio r* versus tau_p for absorber thicknesses 0.75/1.5/2.25 um (mu_n = 500) and electron mobilities
50/100/150 cm2/Vs (W = 1.5 um), with the hole diffusion length at the crossover in units of W.
Input: results/paper_r4_analyses/dimensionless_rstar.csv (cases/analyze_dimensionless_criterion.py). Output: latex/tables/supp_dim_rstar.tex"""
import os
import numpy as np, pandas as pd

root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
T = pd.read_csv(os.path.join(root, "results", "paper_r4_analyses", "dimensionless_rstar.csv"))
cols = [(500, 0.75), (500, 1.5), (500, 2.25), (50, 1.5), (100, 1.5), (150, 1.5)]
taus = [1, 3, 10, 30, 100, 600]
L = [r"\begin{table}[htbp]", r"\centering",
     r"\caption{Crossover ratio $r^*$ (RHJ equal to RBSF--GaP) versus the p-type absorber lifetime for three absorber thicknesses $W$ ($\mu_n=500$~cm$^2$\,V$^{-1}$s$^{-1}$) and three lower electron mobilities ($W=1.5~\mu$m); in parentheses the hole diffusion length of the RHJ absorber at the crossover, $L_h^*/W=\sqrt{\mu_pV_Tr^*\tau_p}/W$ ($\mu_p=112$~cm$^2$\,V$^{-1}$s$^{-1}$). $S_{\rm back}=10^5$~cm\,s$^{-1}$, other parameters nominal; $r^*$ interpolated linearly in $\log r$ (\texttt{cases/analyze\_dimensionless\_criterion.py}).}",
     r"\label{tab:supp-dim}", r"\footnotesize", r"\setlength{\tabcolsep}{3.5pt}",
     r"\begin{tabular}{llrrrrrr}", r"\hline",
     r" & & \multicolumn{3}{c}{$\mu_n=500$, $W$ ($\mu$m)} & \multicolumn{3}{c}{$W=1.5~\mu$m, $\mu_n$} \\",
     r"Condition & $\tau_p$ (ns) & 0.75 & 1.5 & 2.25 & 50 & 100 & 150 \\", r"\hline"]
for cond, cl in (("AM1.5G", "AM1.5G"), ("6500K_1000lx", r"1000~lx")):
    for i, tp in enumerate(taus):
        cells = []
        for mu, W in cols:
            q = T[(T.Condition == cond) & (T.tauP_ns == tp) & (T.muN == mu) & np.isclose(T.W_um, W)]
            cells.append("--" if q.empty or not np.isfinite(q.r_star.iloc[0]) else "%.2f (%.1f)" % (q.r_star.iloc[0], q.Lh_W.iloc[0]))
        L.append("%s & %d & %s \\\\" % (cl if i == 0 else "", tp, " & ".join(cells)))
    if cond == "AM1.5G": L.append(r"\hline")
L += [r"\hline", r"\end{tabular}", r"\end{table}"]
open(os.path.join(root, "latex", "tables", "supp_dim_rstar.tex"), "w", encoding="utf-8", newline=chr(10)).write(chr(10).join(L) + chr(10))
print(chr(10).join(L[10:-3]))
