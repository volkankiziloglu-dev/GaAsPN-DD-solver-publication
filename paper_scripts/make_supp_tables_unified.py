"""Rewrite the rows of tab:supp-edge, tab:supp-edge-eg and tab:supp-thick in latex/supp_mat_en.tex from results/paper_referee_checks_merged
(AM1.5G: earlier runs; 6500 K / 1000 lx: unified extractor). The B_rad column of tab:supp-edge-eg is given in 1e-12 cm3/s. Argument 'write' modifies the file."""
import os, sys
import numpy as np, pandas as pd
root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BS = chr(92)
d = os.path.join(root, "results", "paper_referee_checks_merged")
E = pd.read_csv(os.path.join(d, "referee_checks_E.csv")); T = pd.read_csv(os.path.join(d, "referee_checks_T.csv"))
CB = pd.read_csv(os.path.join(d, "eg_optics_shift_CB.csv")); VB = pd.read_csv(os.path.join(d, "eg_optics_shift_VB.csv"))
print(E.columns.tolist()); print(T.columns.tolist()); print(CB.columns.tolist())
codes = ["HJ-0", "F-HJ", "RBSF-GaP", "RHJ"]


def f(x): return "--" if not np.isfinite(x) else "%.2f" % x


def quad(D, mask_fn, cond):
    return [f(D[mask_fn(D) & (D.Architecture == a) & (D.Condition == cond)].PCE_pct.iloc[0]) for a in codes]


shifts = [(-0.10, "$-100$", "$-$100"), (-0.05, "$-50$", "$-$50"), (0.0, "0", "0"), (0.05, "$+50$", "+50"), (0.10, "$+100$", "+100")]
rows_edge, rows_eg = [], []
for sh, l1, l2 in shifts:
    m = lambda D, sh=sh: abs(D.shift_eV - sh) < 1e-9
    r = quad(E, m, "AM1.5G") + quad(E, m, "6500K_1000lx")
    rows_edge.append("%s & %s %s" % (l1, " & ".join(r), BS + BS))
    r = quad(CB, m, "AM1.5G") + quad(CB, m, "6500K_1000lx")
    brad = CB[m(CB) & (CB.Architecture == "RBSF-GaP") & (CB.Condition == "6500K_1000lx")].BradAbs.iloc[0] / 1e-12
    rows_eg.append("%s & %s & %.2f %s" % (l2, " & ".join(r), brad, BS + BS))
rows_thick = []
for ts, tl in ((1, "$1" + BS + "times$"), (2, "$2" + BS + "times$")):
    for k, (S, sl) in enumerate(((1e5, "$10^{5}$"), (1e3, "$10^{3}$"), (0, "0"))):
        m = lambda D, ts=ts, S=S: (D.thickScale == ts) & (D.Sback == S)
        r = quad(T, m, "AM1.5G") + quad(T, m, "6500K_1000lx")
        rows_thick.append("%s & %s & %s %s" % (tl if k == 0 else "", sl, " & ".join(r), BS + BS))
print("\n".join(rows_edge)); print("\n".join(rows_eg)); print("\n".join(rows_thick))
# max difference CB vs VB
mm = CB.merge(VB, on=["Architecture", "shift_eV", "Condition"], suffixes=("_c", "_v"))
print("max |CB-VB| indoor %.3f, AM1.5G %.3f" % (abs(mm[mm.Condition != "AM1.5G"].PCE_pct_c - mm[mm.Condition != "AM1.5G"].PCE_pct_v).max(),
                                                abs(mm[mm.Condition == "AM1.5G"].PCE_pct_c - mm[mm.Condition == "AM1.5G"].PCE_pct_v).max()))
if len(sys.argv) > 1 and sys.argv[1] == "write":
    fn = os.path.join(root, "latex", "supp_mat_en.tex")
    lines = open(fn, encoding="utf-8").read().split("\n")

    def put(label, rows, first_prefixes):
        i = next(k for k, s in enumerate(lines) if "label{%s}" % label in s)
        j = next(k for k in range(i, len(lines)) if lines[k].startswith(BS + "midrule")) + 1
        for n, rw in enumerate(rows): lines[j + n] = rw
    put("tab:supp-edge", rows_edge, None); put("tab:supp-edge-eg", rows_eg, None); put("tab:supp-thick", rows_thick, None)
    open(fn, "w", encoding="utf-8").write("\n".join(lines)); print("written")
