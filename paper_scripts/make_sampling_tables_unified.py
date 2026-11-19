"""Regenerate the sampling tables of the manuscript and the Supplementary Material from the merged data (results/*_merged; AM1.5G rows
from the *_unified_am reruns, indoor rows from the *_unified reruns, both with the unified extractor; see cases/merge_unified.py).
Valid sample: AM1.5G = class AM15G and finite PCE; indoor = finite PCE (the acceptance class is a label only).
Percentiles: NumPy linear interpolation (type 7) for every table (2026-11-07; the nanosecond-sampling table used the MATLAB/Hazen
estimator mprc() before, kept below for reference only).
Usage: py make_sampling_tables_unified.py [write]"""
import os, sys
import numpy as np, pandas as pd
from scipy.stats import spearmanr
root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = os.path.join(root, "results"); BS = chr(92)
A = ["HJ-0", "F-HJ", "RBSF-GaP", "RHJ"]
S1 = pd.read_csv(os.path.join(R, "paper_uncertainty_merged", "uncertainty_samples.csv"))
S2 = pd.read_csv(os.path.join(R, "paper_uncertainty_ns_merged", "ns_samples_all.csv"))
BT = pd.read_csv(os.path.join(R, "paper_uncertainty_merged", "spearman_bootstrap.csv"))
WL = pd.read_csv(os.path.join(R, "paper_uncertainty_merged", "win_tie_loss.csv"))
conds = ["AM1.5G", "6500K_1000lx"]


def valid(t, cond):
    ok = np.isfinite(t.PCE_pct)
    return ok & (t.Class == "AM15G") if cond == "AM1.5G" else ok


def pivot(T, cond):
    t = T[T.Condition == cond]; t = t[valid(t, cond)]
    return t.pivot(index="Sample", columns="Architecture", values="PCE_pct").reindex(range(1, 61))


def sg(x, nd=2):
    s = ("%." + str(nd) + "f") % abs(x)
    return ("$-$" + s) if (x < 0 and float(s) != 0) else s


def sgm(x, nd=2):   # math-mode minus, as used inside tables of the supplement
    s = ("%." + str(nd) + "f") % abs(x)
    return ("$-" + s + "$") if (x < 0 and float(s) != 0) else s


def mprc(v, p):
    s = np.sort(np.asarray(v, float)); n = len(s); r = []
    for pp in p:
        x = pp / 100 * n + 0.5
        if x <= 1: r.append(s[0])
        elif x >= n: r.append(s[-1])
        else:
            k = int(np.floor(x)); r.append(s[k - 1] + (x - k) * (s[k] - s[k - 1]))
    return r


def f2(x): 
    s = "%.2f" % x
    return "0.00" if float(s) == 0 else s


out = {}
# ---- main table tab:uncertainty-pce
nomAM = {"HJ-0": 19.34, "F-HJ": 21.83, "RBSF-GaP": 22.36, "RHJ": 22.21}
nomIN = {"HJ-0": 40.66, "F-HJ": 45.03, "RBSF-GaP": 45.45, "RHJ": 45.15}
nm = {"HJ-0": "HJ-0", "F-HJ": "F-HJ", "RBSF-GaP": "RBSF--GaP", "RHJ": "RHJ"}
rows = []
P1 = {c: pivot(S1, c) for c in conds}
for a in A:
    cells = []
    for c, nom in zip(conds, (nomAM, nomIN)):
        v = P1[c][a].dropna(); q = np.percentile(v, [5, 50, 95])
        cells += ["%d" % len(v), "%.2f" % q[0], "%.2f" % q[1], "%.2f" % q[2], "%.2f" % nom[a]]
    rows.append("%s & %s %s" % (nm[a], " & ".join(cells), BS + BS))
out["main1"] = rows
# ---- supp: ranking frequencies sampling 1
pairs = [("RBSF-GaP", "RHJ"), ("RBSF-GaP", "HJ-0"), ("RHJ", "HJ-0"), ("F-HJ", "HJ-0"), ("RBSF-GaP", "F-HJ"), ("RHJ", "F-HJ")]
rows = []
for a, b in pairs:
    cells = []
    for c in conds:
        d = (P1[c][a] - P1[c][b]).dropna(); q = np.percentile(d, [5, 50, 95])
        cells += ["%d/%d" % ((d > 0).sum(), len(d)), "/".join(sgm(x) for x in q)]
    rows.append("%s $>$ %s & %s %s" % (nm[a], nm[b], " & ".join(cells), BS + BS))
out["order1"] = rows
# ---- supp: Spearman sampling 1 (samples valid for all four)
P = ["tauAbs", "Sint", "Nabs", "chiGaP", "tauGaP", "muN"]
rows = []
for ci, c in enumerate(conds):
    piv = P1[c]; full = piv.dropna(subset=A); prm = S1[S1.Condition == c].drop_duplicates("Sample").set_index("Sample").loc[full.index]
    lab = ["PCE, HJ-0", "PCE, F-HJ", "PCE, RBSF--GaP", "PCE, RHJ", "RBSF--GaP$-$RHJ"]
    for k, nmq in enumerate(lab):
        y = full[A[k]] if k < 4 else full["RBSF-GaP"] - full["RHJ"]
        rr = [spearmanr(prm[p], y)[0] for p in P]
        first = ("AM1.5G" if ci == 0 else "6500~K, 1000~lx") if k == 0 else ""
        rows.append("%s & %s & %s %s" % (first, nmq, " & ".join(sgm(x) for x in rr), BS + BS))
    if ci == 0: rows.append(BS + "midrule")
out["rho1"] = rows
# ---- supp: sampling 2 PCE distributions and Spearman
P2 = {c: pivot(S2, c) for c in conds}
rows = []
for a in A:
    cells = []
    for c in conds:
        v = P2[c][a].dropna(); q = np.percentile(v, [5, 50, 95])
        cells += ["%d" % len(v), "%.2f/%.2f/%.2f" % tuple(q), ""]
    rows.append("%s & %s %s" % (nm[a], " & ".join(cells), BS + BS))
out["ns_pce"] = rows
P8 = ["tauP", "ratio", "Sback", "Sint", "Nabs", "chiGaP", "tauGaP", "muN"]
rows = []
for ci, c in enumerate(conds):
    piv = P2[c]; prm = S2[S2.Condition == c].drop_duplicates("Sample").set_index("Sample").reindex(range(1, 61))
    lab = ["PCE, HJ-0", "PCE, F-HJ", "PCE, RBSF--GaP", "PCE, RHJ", "RBSF--GaP$-$RHJ"]
    for k, nmq in enumerate(lab):
        if k < 4: y = piv[A[k]]
        else: y = piv["RBSF-GaP"] - piv["RHJ"]
        m = y.notna()
        rr = [spearmanr(prm.loc[m, p], y[m])[0] for p in P8]
        first = ("AM1.5G" if ci == 0 else "6500~K, 1000~lx") if k == 0 else ""
        rows.append("%s & %s & %s %s" % (first, nmq, " & ".join(sgm(x) for x in rr), BS + BS))
    if ci == 0: rows.append(BS + "midrule")
out["rho2"] = rows
# ---- bootstrap table
spec = [("set1", "PCE F-HJ vs tauAbs", "PCE, F-HJ, vs $" + BS + "tau_{" + BS + "mathrm{abs}}$"), ("set1", "PCE RBSF-GaP vs tauAbs", "PCE, RBSF--GaP, vs $" + BS + "tau_{" + BS + "mathrm{abs}}$"),
        ("set1", "PCE RHJ vs tauAbs", "PCE, RHJ, vs $" + BS + "tau_{" + BS + "mathrm{abs}}$"), ("set1", "PCE HJ-0 vs Nabs", "PCE, HJ-0, vs $N_{" + BS + "mathrm{abs}}$"),
        ("set1", "RBSF-RHJ vs tauAbs", "RBSF--GaP$-$RHJ vs $" + BS + "tau_{" + BS + "mathrm{abs}}$"),
        ("set2", "RBSF-RHJ vs ratio", "RBSF--GaP$-$RHJ vs $r$"), ("set2", "RBSF-RHJ vs tauP", "RBSF--GaP$-$RHJ vs $" + BS + "tau_p$"), ("set2", "PCE RHJ vs ratio", "PCE, RHJ, vs $r$"),
        ("set2", "PCE RBSF-GaP vs tauP", "PCE, RBSF--GaP, vs $" + BS + "tau_p$"), ("set2", "PCE HJ-0 vs Sback", "PCE, HJ-0, vs $S_{" + BS + "mathrm{back}}$"), ("set2", "RBSF-RHJ vs Sback", "RBSF--GaP$-$RHJ vs $S_{" + BS + "mathrm{back}}$")]
rows = []
for k, (st, q, lab) in enumerate(spec):
    cells = []
    for c in conds:
        r = BT[(BT.set == st) & (BT.condition == c) & (BT.quantity == q)].iloc[0]
        cells.append("%s [%s, %s] (%d)" % (sg(r.rho), sg(r.lo95), sg(r.hi95), r.n))
    rows.append("%d & %s & %s %s" % (1 if st == "set1" else 2, lab, " & ".join(cells), BS + BS))
    if k == 4: rows.append(BS + "midrule")
out["boot"] = rows
# ---- win/tie/loss table
rows = []
for a, b in [("RBSF-GaP", "RHJ"), ("RBSF-GaP", "F-HJ"), ("RBSF-GaP", "HJ-0"), ("RHJ", "F-HJ"), ("RHJ", "HJ-0")]:
    for ci, c in enumerate(conds):
        cells = []
        for st in (1, 2):
            r = WL[(WL.set == st) & (WL.condition == c) & (WL.A == a) & (WL.B == b)].iloc[0]
            cells += ["%d/%d/%d" % (r.A_wins, r.ties, r.A_loses), "%s [%s, %s]" % (f2(r.p_win), f2(max(r.win_lo, 0)), f2(r.win_hi))]
        first = ("%s $>$ %s" % (nm[a], nm[b])) if ci == 0 else ""
        rows.append("%s & %s & %s %s" % (first, "AM1.5G" if ci == 0 else "1000~lx", " & ".join(cells), BS + BS))
out["wtl"] = rows
for k, v in out.items():
    print("== " + k); print("\n".join(v))

if len(sys.argv) > 1 and sys.argv[1] == "write":
    def block(lines, start_pat, end_pat, new_rows, from_start=False):
        i = next(k for k, s in enumerate(lines) if start_pat in s)
        j = next(k for k in range(i, len(lines)) if end_pat in lines[k])
        lines[i:j] = new_rows
    mf = os.path.join(root, "latex", "manuscript_en.tex"); sf = os.path.join(root, "latex", "supp_mat_en.tex")
    ml = open(mf, encoding="utf-8", newline="").read().split("\n"); sl = open(sf, encoding="utf-8", newline="").read().split("\n")   # keep LF line endings
    block(ml, "HJ-0 & 60 & 18.89", BS + "bottomrule", out["main1"])
    # supplementary blocks: replace from the first data row to the line before bottomrule of the table that follows the label
    def repl_after(label, first_row_start, rows):
        i = next(k for k, s in enumerate(sl) if "label{%s}" % label in s)
        j = next(k for k in range(i, len(sl)) if sl[k].startswith(first_row_start))
        e = next(k for k in range(j, len(sl)) if sl[k].startswith(BS + "bottomrule"))
        sl[j:e] = rows
    repl_after("tab:supp-uncertainty-order", "RBSF--GaP $>$ RHJ &", out["order1"])
    repl_after("tab:supp-uncertainty-rho", "AM1.5G & PCE, HJ-0", out["rho1"])
    repl_after("tab:supp-ns-pce", "HJ-0 &", out["ns_pce"])
    repl_after("tab:supp-ns-rho", "AM1.5G & PCE, HJ-0", out["rho2"])
    repl_after("tab:supp-bootstrap", "1 & PCE, F-HJ", out["boot"])
    repl_after("tab:supp-wtl", "RBSF--GaP $>$ RHJ &", out["wtl"])
    open(mf, "w", encoding="utf-8", newline="").write("\n".join(ml)); open(sf, "w", encoding="utf-8", newline="").write("\n".join(sl)); print("written")
