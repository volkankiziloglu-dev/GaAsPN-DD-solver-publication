"""Win/tie/loss frequencies with an equivalence margin, Wilson intervals, imputation bounds and a failure-versus-parameter table for the two 60-sample parameter samplings.
Valid sample per architecture: AM1.5G = class AM15G and finite PCE; 6500 K / 1000 lx = finite PCE (unified extractor; the acceptance class is a label only).
Margin: |difference| < 0.1 percentage points is a tie (several times the voltage-step error of the unified extractor, 0.014 points AM1.5G and 0.021 points at 1000 lx; sensitivity to 0.05 and 0.2 points in analyze_sampling_r4_minor.py).
Outputs: results/paper_uncertainty_merged/win_tie_loss.csv, results/paper_uncertainty_merged/failure_vs_parameters.csv
"""
import os
import numpy as np, pandas as pd
from scipy.stats import beta, mannwhitneyu

root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
A = ["HJ-0", "F-HJ", "RBSF-GaP", "RHJ"]
MARGIN = 0.1


def wilson(k, n, z=1.96):
    if n == 0: return (np.nan, np.nan)
    p = k / n; d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d; h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return c - h, c + h


def load(setno):
    if setno == 1:
        T = pd.read_csv(os.path.join(root, "results", "paper_uncertainty_merged", "uncertainty_samples.csv"))
        P = ["tauAbs", "Sint", "Nabs", "chiGaP", "tauGaP", "muN"]
    else:
        T = pd.read_csv(os.path.join(root, "results", "paper_uncertainty_ns_merged", "ns_samples_all.csv"))
        P = ["tauP", "ratio", "Sback", "Sint", "Nabs", "chiGaP", "tauGaP", "muN"]
    return T, P


def valid_mask(t, cond, setno, allow_screen):
    fin = np.isfinite(t.PCE_pct)
    if cond == "AM1.5G":
        return fin & (t.Class == "AM15G")
    return fin   # indoor (unified extractor): all finite points; the acceptance class is a label only


rows = []; frows = []
pairs = [("RBSF-GaP", "RHJ"), ("RBSF-GaP", "F-HJ"), ("RBSF-GaP", "HJ-0"), ("RHJ", "F-HJ"), ("RHJ", "HJ-0"), ("F-HJ", "HJ-0")]
for setno in (1, 2):
    T, P = load(setno)
    nS = T.Sample.nunique()
    for cond in ["AM1.5G", "6500K_1000lx"]:
        t = T[T.Condition == cond]
        ok = valid_mask(t, cond, setno, allow_screen=(setno == 1))
        piv = t[ok].pivot(index="Sample", columns="Architecture", values="PCE_pct").reindex(range(1, nS + 1))
        for a, b in pairs:
            both = piv[[a, b]].dropna()
            d = both[a] - both[b]; n = len(d)
            w = int((d >= MARGIN).sum()); l = int((d <= -MARGIN).sum()); tie = n - w - l
            lo, hi = wilson(w, n); lo_l, hi_l = wilson(l, n)
            miss = nS - n
            rows.append((setno, cond, a, b, n, miss, w, tie, l, w / n, lo, hi, l / n, lo_l, hi_l, w / nS, (w + miss) / nS))
        # failure versus parameters (a sample is incomplete for the pair RBSF-GaP/RHJ if either is invalid)
        prm = T.drop_duplicates("Sample").set_index("Sample")[P].reindex(range(1, nS + 1))
        bad = piv[["RBSF-GaP", "RHJ"]].isna().any(axis=1)
        for p in P:
            x = prm[p]; xb, xg = x[bad], x[~bad]
            if len(xb) >= 3 and len(xg) >= 3:
                u = mannwhitneyu(xb, xg).pvalue
                # fraction incomplete in top and bottom tercile of the parameter (log for positive-valued)
                q = np.log10(x) if p not in ("chiGaP",) else x
                tercile = pd.qcut(q, 3, labels=["low", "mid", "high"])
                fr = bad.groupby(tercile, observed=True).mean()
                frows.append((setno, cond, p, int(bad.sum()), nS, fr["low"], fr["mid"], fr["high"], u))
cols = ["set", "condition", "A", "B", "n_valid", "n_missing", "A_wins", "ties", "A_loses", "p_win", "win_lo", "win_hi", "p_loss", "loss_lo", "loss_hi", "win_lower_bound", "win_upper_bound"]
R = pd.DataFrame(rows, columns=cols)
R.to_csv(os.path.join(root, "results", "paper_uncertainty_merged", "win_tie_loss.csv"), index=False, float_format="%.3f")
F = pd.DataFrame(frows, columns=["set", "condition", "parameter", "n_incomplete", "n_samples", "frac_low", "frac_mid", "frac_high", "MannWhitney_p"])
F.to_csv(os.path.join(root, "results", "paper_uncertainty_merged", "failure_vs_parameters.csv"), index=False, float_format="%.3f")
pd.set_option("display.width", 220)
print(R.to_string(index=False, float_format=lambda v: "%.2f" % v))
print(F.to_string(index=False, float_format=lambda v: "%.3f" % v))
