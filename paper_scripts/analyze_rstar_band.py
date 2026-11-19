"""Indifference band of the crossover ratio: the ratios r at which PCE_RHJ - PCE_RBSF = -m and +m
(m = 0.1 points, the equivalence margin; 0.2 points in the nanosecond regime), linear interpolation in log r on the design-map grid.
Input: results/paper_design_map_merged/design_map_[abcd].csv (S_back = 1e5, mu_n = 500, W = 1.5 um). Output: results/paper_design_map_merged/rstar_band.csv"""
import glob, os
import numpy as np, pandas as pd
root = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); R = os.path.join(root, "results", "paper_design_map_merged")
D = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(os.path.join(R, "design_map_[abcd].csv")))]); D = D[D.Sback == 1e5]


def cross(r, d, level):
    lr = np.log10(r); k = np.where(d >= level)[0]
    if len(k) == 0: return np.inf
    k = k[0]
    if k == 0: return r[0] if d[0] == level else -r[0]
    return 10 ** (lr[k - 1] + (level - d[k - 1]) * (lr[k] - lr[k - 1]) / (d[k] - d[k - 1]))


rows = []
for (cond, tp), g in D.groupby(["Condition", "tauP_ns"]):
    rb = g[g.Architecture == "RBSF-GaP"].PCE_pct.iloc[0]
    rh = g[g.Architecture == "RHJ"].sort_values("ratio"); r = rh.ratio.values; d = rh.PCE_pct.values - rb
    rows.append((cond, tp, cross(r, d, -0.1), cross(r, d, 0.0), cross(r, d, 0.1), cross(r, d, -0.2), cross(r, d, 0.2)))
T = pd.DataFrame(rows, columns=["Condition", "tauP_ns", "r_m01", "r_star", "r_p01", "r_m02", "r_p02"])
T.to_csv(os.path.join(R, "rstar_band.csv"), index=False, lineterminator=chr(10), float_format="%.3f")
print(T.round(2).to_string(index=False))
