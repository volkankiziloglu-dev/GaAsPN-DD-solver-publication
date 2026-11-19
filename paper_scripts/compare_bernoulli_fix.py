"""Compare the BERN re-runs (sgBernoulli with expm1, 2026-11-08) with the original OPT and L200 rows (x/(exp(x)-1)).
Output: results/paper_r4_analyses/bernoulli_fix_comparison.csv; prints max |dPCE| per condition and class changes."""
import glob, os
import pandas as pd

root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = os.path.join(root, "results", "paper_r4_analyses")
rd = lambda pat: pd.concat([pd.read_csv(f) for f in sorted(glob.glob(os.path.join(R, pat)))], ignore_index=True)
old = pd.concat([rd("OPT_*of4.csv"), rd("L200_*of2.csv")], ignore_index=True)
new = rd("BERN_*of4.csv")
key = ["Architecture", "tauP_ns", "ratio", "Sback", "tScale", "Nabs", "Condition"]
for d in (old, new):
    for k in ["tauP_ns", "ratio", "Sback", "tScale", "Nabs"]:
        d[k] = d[k].astype(float)
m = new.merge(old, on=key, suffixes=("_new", "_old"), how="left")
for q in ["Jsc", "Voc_mV", "FF_pct", "PCE_pct"]:
    m["d_" + q] = m[q + "_new"] - m[q + "_old"]
cols = key + ["PCE_pct_old", "PCE_pct_new", "d_PCE_pct", "d_Jsc", "d_Voc_mV", "d_FF_pct", "Class_old", "Class_new"]
m = m[cols]
m.to_csv(os.path.join(R, "bernoulli_fix_comparison.csv"), index=False, lineterminator=chr(10))
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 100)
print(m.to_string(index=False))
print()
print(m.groupby("Condition")["d_PCE_pct"].agg(lambda s: s.abs().max()))
print("class changes:", int((m["Class_old"] != m["Class_new"]).sum()), "of", len(m))
