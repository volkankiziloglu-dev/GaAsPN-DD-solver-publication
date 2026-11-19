"""Supplementary table: RBSF-GaP - RHJ (PCE points) with architecture-specific leakage at 100-1000 lx, 0.6 us and 10 ns.
Input: results/paper_leakage_v112/leakage_archspecific_v112_wide.csv (paper_scripts/analyze_leakage_v112.py),
       rows Set == "voc_lux" and HJbranch == "old" (homojunction branch J0p = 4.97e-10 A/cm2, n = 3.72, Rp = 8.82e3 Ohm cm2).
Output: latex/tables/supp_leak_archspecific_lux.tex"""
import os
import numpy as np
import pandas as pd

root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
T = pd.read_csv(os.path.join(root, "results", "paper_leakage_v112", "leakage_archspecific_v112_wide.csv"))
T = T[(T.Set == "voc_lux") & (T.HJbranch == "old")]
L = [r"\begin{table}[htbp]", r"\centering",
     r"\caption{RBSF--GaP$-$RHJ (PCE points, CIE 6500~K) with architecture-specific leakage: the branch fitted to the Klitzke homojunction for HJ-0, F-HJ and RBSF--GaP and the measured excess current of the Klitzke RHJ (relative to the model with $N_D=2\times10^{16}$~cm$^{-3}$, digitized at 100, 200, 500 and 1000~lx) for the RHJ, both scaled by the factor given. Leading design in parentheses where it is not RBSF--GaP. Post-processing of the intrinsic curves.}",
     r"\label{tab:supp-leaklux}", r"\small",
     r"\begin{tabular}{rrrrrr}", r"\hline",
     r"$\tau_p$ & Illuminance & Intrinsic & $\times0.5$ & $\times1$ & $\times2$ \\", r"\hline"]


def fmt(d):
    d = float(np.round(d, 2))
    s = f"{abs(d):.2f}"
    return ("+" if d >= 0 else r"$-$") + s


for tau, tl in ((600, r"0.6~$\mu$s"), (10, r"10~ns")):
    for i, lux in enumerate((1000, 500, 200, 100)):
        cells = []
        for sc in (0, 0.5, 1, 2):
            r = T[(T.tauP_ns == tau) & (T.Lux == lux) & ((T.Scale - sc).abs() < 1e-9)].iloc[0]
            c = fmt(r.RBSF_minus_RHJ)
            if r.best != "RBSF-GaP":
                c += f" ({r.best})"
            cells.append(c)
        L.append(f"{tl if i == 0 else ''} & {lux}~lx & " + " & ".join(cells) + r" \\")
    L.append(r"\hline")
L += [r"\end{tabular}", r"\end{table}"]
out = os.path.join(root, "latex", "tables", "supp_leak_archspecific_lux.tex")
os.makedirs(os.path.dirname(out), exist_ok=True)
open(out, "w", encoding="utf-8", newline="\n").write("\n".join(L) + "\n")
print(out)
