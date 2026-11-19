"""Small statistical checks requested by referee R4 (third simulated review, docs/HAKEM_SIMULASYONU_2026-11-06.md).
1. Realized input-input Spearman correlations of the two Latin hypercube samplings (no correlation control was used).
2. Partial rank correlations (PRCC: Spearman correlation of the rank residuals after regressing out the other sampled inputs) for the
   coefficients quoted in the text, to show whether input correlations confound them.
3. Win/tie/loss counts with equivalence margins of 0.05, 0.1 (used) and 0.2 points.
4. Percentiles with one estimator for all tables: NumPy linear interpolation (type 7), the estimator of Table tab:uncertainty-pce.
Inputs : results/paper_uncertainty_merged/uncertainty_samples.csv, results/paper_uncertainty_ns_merged/ns_samples_all.csv
Output : results/paper_uncertainty_merged/r4_minor_checks.txt (printed) and r4_input_correlations.csv, r4_prcc.csv, r4_wtl_margins.csv
"""
import os
import numpy as np, pandas as pd
from scipy.stats import spearmanr, rankdata

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = os.path.join(ROOT, "results")
OUT = os.path.join(R, "paper_uncertainty_merged")
S1 = pd.read_csv(os.path.join(R, "paper_uncertainty_merged", "uncertainty_samples.csv"))
S2 = pd.read_csv(os.path.join(R, "paper_uncertainty_ns_merged", "ns_samples_all.csv"))
P1 = ["tauAbs", "Sint", "Nabs", "chiGaP", "tauGaP", "muN"]
P2 = ["tauP", "ratio", "Sback", "Sint", "Nabs", "chiGaP", "tauGaP", "muN"]
CONDS = ["AM1.5G", "6500K_1000lx"]
lines = []
log = lambda s="": (print(s), lines.append(s))


def valid(t, cond):
    ok = np.isfinite(t.PCE_pct)
    return ok & (t.Class == "AM15G") if cond == "AM1.5G" else ok


def pivot(T, cond):
    t = T[(T.Condition == cond)]; t = t[valid(t, cond)]
    return t.pivot(index="Sample", columns="Architecture", values="PCE_pct").reindex(range(1, 61))


def params(T, P):
    return T.drop_duplicates("Sample").set_index("Sample").reindex(range(1, 61))[P]


# ---- 1. input correlations
rows = []
for name, T, P in (("sampling 1", S1, P1), ("sampling 2", S2, P2)):
    X = params(T, P)
    rho = spearmanr(X.values).correlation
    off = [(P[i], P[j], rho[i, j]) for i in range(len(P)) for j in range(i + 1, len(P))]
    off.sort(key=lambda z: -abs(z[2]))
    log("%s: max |rho| between inputs = %.3f (%s-%s); next %s" % (name, abs(off[0][2]), off[0][0], off[0][1],
        ", ".join("%s-%s %.2f" % (a, b, r) for a, b, r in off[1:4])))
    rows += [dict(sampling=name, a=a, b=b, rho=round(r, 4)) for a, b, r in off]
pd.DataFrame(rows).to_csv(os.path.join(OUT, "r4_input_correlations.csv"), index=False)


# ---- 2. PRCC
def prcc(y, X, j):
    Ry = rankdata(y); RX = np.column_stack([rankdata(X[:, k]) for k in range(X.shape[1])])
    others = np.column_stack([np.ones(len(y))] + [RX[:, k] for k in range(X.shape[1]) if k != j])
    ry = Ry - others @ np.linalg.lstsq(others, Ry, rcond=None)[0]
    rx = RX[:, j] - others @ np.linalg.lstsq(others, RX[:, j], rcond=None)[0]
    return np.corrcoef(ry, rx)[0, 1]


wanted = [("sampling 1", S1, P1, "RBSF-GaP", "RHJ", "tauAbs"), ("sampling 1", S1, P1, "RBSF-GaP", "RHJ", "chiGaP"),
          ("sampling 1", S1, P1, "HJ-0", None, "Nabs"), ("sampling 1", S1, P1, "RBSF-GaP", None, "Nabs"),
          ("sampling 2", S2, P2, "RBSF-GaP", "RHJ", "ratio"), ("sampling 2", S2, P2, "RBSF-GaP", "RHJ", "tauP"),
          ("sampling 2", S2, P2, "RBSF-GaP", "RHJ", "Sback"), ("sampling 2", S2, P2, "HJ-0", None, "Sback"),
          ("sampling 2", S2, P2, "HJ-0", None, "Sint"), ("sampling 2", S2, P2, "RHJ", None, "ratio")]
prow = []
for name, T, P, a, b, par in wanted:
    X = params(T, P)
    for c in CONDS:
        piv = pivot(T, c)
        y = piv[a] - piv[b] if b else piv[a]
        ok = np.isfinite(y.values)
        rs = spearmanr(y.values[ok], X[par].values[ok]).correlation
        pr = prcc(y.values[ok], X.values[ok], P.index(par))
        lab = ("%s-%s" % (a, b)) if b else ("PCE %s" % a)
        prow.append(dict(sampling=name, condition=c, quantity=lab, parameter=par, spearman=round(rs, 3), prcc=round(pr, 3), n=int(ok.sum())))
        log("%s %-13s %-18s vs %-6s Spearman %+.3f  PRCC %+.3f  (n=%d)" % (name, c, lab, par, rs, pr, ok.sum()))
pd.DataFrame(prow).to_csv(os.path.join(OUT, "r4_prcc.csv"), index=False)

# ---- 3. win/tie/loss vs margin
wrow = []
for name, T in (("sampling 1", S1), ("sampling 2", S2)):
    for c in CONDS:
        piv = pivot(T, c)
        for a, b in (("RBSF-GaP", "RHJ"), ("RBSF-GaP", "F-HJ"), ("RBSF-GaP", "HJ-0"), ("RHJ", "F-HJ")):
            d = (piv[a] - piv[b]).dropna().values
            for m in (0.05, 0.1, 0.2):
                w, t, l = int((d >= m).sum()), int((abs(d) < m).sum()), int((d <= -m).sum())
                wrow.append(dict(sampling=name, condition=c, A=a, B=b, margin=m, wins=w, ties=t, losses=l, n=len(d)))
            log("%s %-13s %-9s vs %-5s  W/T/L  m=0.05: %s  m=0.1: %s  m=0.2: %s" % (name, c, a, b,
                *["%d/%d/%d" % (r["wins"], r["ties"], r["losses"]) for r in wrow[-3:]]))
pd.DataFrame(wrow).to_csv(os.path.join(OUT, "r4_wtl_margins.csv"), index=False)

# ---- 4. percentiles (NumPy linear) for the nanosecond tables
log()
for c in CONDS:
    piv = pivot(S2, c)
    for a in ["HJ-0", "F-HJ", "RBSF-GaP", "RHJ"]:
        v = piv[a].dropna().values
        log("ns %-13s PCE %-9s P5/P50/P95 linear %.2f/%.2f/%.2f (n=%d)" % (c, a, *np.percentile(v, [5, 50, 95]), len(v)))
    for a, b in (("RBSF-GaP", "RHJ"), ("RBSF-GaP", "F-HJ"), ("RBSF-GaP", "HJ-0"), ("F-HJ", "HJ-0"), ("RHJ", "HJ-0")):
        d = (piv[a] - piv[b]).dropna().values
        log("ns %-13s %s>%s  %d/%d  median (P5, P95) linear %.2f (%.2f, %.2f)" % (c, a, b, (d > 0).sum(), len(d), *np.percentile(d, [50, 5, 95])))

n = 60; z = 1.96
for p in (0.5,):
    hw = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    log("Wilson half-width at p=0.5, n=60: %.3f" % hw)
open(os.path.join(OUT, "r4_minor_checks.txt"), "w").write("\n".join(lines) + "\n")
