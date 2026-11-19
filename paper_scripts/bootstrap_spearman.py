"""Bootstrap (percentile, 95 %) confidence intervals for the Spearman correlations quoted in the manuscript.
Set 1: results/paper_uncertainty_merged/uncertainty_samples.csv  (60 samples; AM1.5G class AM15G, indoor all finite PCE (unified extractor; the class is a label only); correlations on samples valid for all four architectures, as in analyze_uncertainty_sampling.js)
Set 2: results/paper_uncertainty_ns_merged/ns_samples_all.csv   (60 samples, ns regime; correlations per architecture on its valid samples, RBSF-GaP minus RHJ on samples valid for both, as in analyze_uncertainty_sampling_ns.m)
Resampling is over samples (rows of the sample table), 10000 resamples, seed 20261105. Output: results/paper_uncertainty_merged/spearman_bootstrap.csv
"""
import os, sys
import numpy as np, pandas as pd
from scipy.stats import spearmanr

root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
rng = np.random.default_rng(20261105)
NB = 10000
A = ["HJ-0", "F-HJ", "RBSF-GaP", "RHJ"]


def boot(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float); n = len(x)
    r0 = spearmanr(x, y)[0]
    idx = rng.integers(0, n, size=(NB, n))
    rs = np.empty(NB)
    for b in range(NB):
        i = idx[b]
        rs[b] = spearmanr(x[i], y[i])[0]
    rs = rs[np.isfinite(rs)]
    return r0, np.percentile(rs, 2.5), np.percentile(rs, 97.5), n


rows = []

# ---------------- set 1
T = pd.read_csv(os.path.join(root, "results", "paper_uncertainty_merged", "uncertainty_samples.csv"))
for cond in ["AM1.5G", "6500K_1000lx"]:
    t = T[T.Condition == cond]
    if cond == "AM1.5G":
        ok = t.Class == "AM15G"
    else:
        ok = np.isfinite(t.PCE_pct)   # indoor (unified extractor): all finite points, class is a label only
    t = t[ok & np.isfinite(t.PCE_pct)]
    piv = t.pivot(index="Sample", columns="Architecture", values="PCE_pct")
    full = piv.dropna(subset=A)
    prm = t.drop_duplicates("Sample").set_index("Sample").loc[full.index]
    d = full["RBSF-GaP"] - full["RHJ"]
    for a in A:
        for p in ["tauAbs", "Nabs"]:
            r = boot(prm[p], full[a]); rows.append(("set1", cond, f"PCE {a} vs {p}") + r)
    r = boot(prm["tauAbs"], d); rows.append(("set1", cond, "RBSF-RHJ vs tauAbs") + r)

# ---------------- set 2
T = pd.read_csv(os.path.join(root, "results", "paper_uncertainty_ns_merged", "ns_samples_all.csv"))
for cond in ["AM1.5G", "6500K_1000lx"]:
    t = T[T.Condition == cond]
    ok = np.isfinite(t.PCE_pct)
    if cond != "AM1.5G":
        pass  # unified extractor: indoor uses all finite points
    t = t[ok]
    piv = t.pivot(index="Sample", columns="Architecture", values="PCE_pct")
    prm = T.drop_duplicates("Sample").set_index("Sample")
    both = piv.dropna(subset=["RBSF-GaP", "RHJ"])
    d = both["RBSF-GaP"] - both["RHJ"]
    for p in ["ratio", "tauP", "Sback"]:
        r = boot(prm.loc[both.index, p], d); rows.append(("set2", cond, f"RBSF-RHJ vs {p}") + r)
    for a in A:
        v = piv[a].dropna()
        r = boot(prm.loc[v.index, "tauP"], v); rows.append(("set2", cond, f"PCE {a} vs tauP") + r)
    v = piv["RHJ"].dropna(); r = boot(prm.loc[v.index, "ratio"], v); rows.append(("set2", cond, "PCE RHJ vs ratio") + r)
    v = piv["HJ-0"].dropna(); r = boot(prm.loc[v.index, "Sback"], v); rows.append(("set2", cond, "PCE HJ-0 vs Sback") + r)

out = pd.DataFrame(rows, columns=["set", "condition", "quantity", "rho", "lo95", "hi95", "n"])
out.to_csv(os.path.join(root, "results", "paper_uncertainty_merged", "spearman_bootstrap.csv"), index=False, float_format="%.3f")
pd.set_option("display.width", 200)
print(out.to_string(index=False, float_format=lambda v: "%.3f" % v))
