"""Supplementary table: RBSF-GaP - RHJ (PCE points) with architecture-specific leakage at 100-1000 lx, 0.6 us and 10 ns.
Input: results/paper_leakage_robustness/leakage_archspecific_lux.csv (cases/analyze_leakage_archspecific_lux.py).
Output: latex/tables/supp_leak_archspecific_lux.tex"""
import os
import pandas as pd

root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
T = pd.read_csv(os.path.join(root, "results", "paper_leakage_robustness", "leakage_archspecific_lux.csv"))
T = T[T.Mode == "archspecific"]
L = [r"\begin{table}[htbp]", r"\centering",
     r"\caption{RBSF--GaP$-$RHJ (PCE points, CIE 6500~K) with architecture-specific leakage: the branch fitted to the Klitzke homojunction for HJ-0, F-HJ and RBSF--GaP and the measured excess current of the Klitzke RHJ for the RHJ, both scaled by the factor given; at 200 and 500~lx the RHJ excess current is interpolated per voltage as a power law in illuminance between the measurements at 100 and 1000~lx. Leading design in parentheses where it is not RBSF--GaP. Post-processing of the intrinsic curves (\texttt{cases/analyze\_leakage\_archspecific\_lux.py}).}",
     r"\label{tab:supp-leaklux}", r"\small",
     r"\begin{tabular}{rrrrrr}", r"\hline",
     r"$\tau_p$ & Illuminance & Intrinsic & $\times0.5$ & $\times1$ & $\times2$ \\", r"\hline"]
for tau, tl in ((600, r"0.6~$\mu$s"), (10, r"10~ns")):
    for i, lux in enumerate((1000, 500, 200, 100)):
        cells = []
        for sc in (0, 0.5, 1, 2):
            g = T[(T.tauP_ns == tau) & (T.Lux == lux) & (T.Scale == sc)].set_index("Architecture").PCE_pct
            d = g["RBSF-GaP"] - g["RHJ"]; best = g.idxmax()
            cells.append(("%+.2f" % d).replace("-", "$-$") +("" if best == "RBSF-GaP" else " (%s)" % best))
        L.append("%s & %d~lx & %s \\\\" % (tl if i == 0 else "", lux, " & ".join(cells)))
    if tau == 600: L.append(r"\hline")
L += [r"\hline", r"\end{tabular}", r"\end{table}"]
open(os.path.join(root, "latex", "tables", "supp_leak_archspecific_lux.tex"), "w", encoding="utf-8", newline=chr(10)).write(chr(10).join(L) + chr(10))
print(chr(10).join(L[8:-3]))
