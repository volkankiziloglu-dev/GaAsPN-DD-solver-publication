"""Prior sensitivity of the nanosecond-regime ranking frequency.

The RBSF-GaP vs RHJ sign is set almost entirely by r = tau_n/tau_p relative to the deterministic crossover r*(tau_p) of the design map
(RHJ ahead iff r > r*). Hence the sample frequency "RHJ ahead" of the ns sampling is P(r > r*(tau_p)) under the priors.
1. Check: classify each ns-sampling pair by r > r*(tau_p) (log-log interpolation of the design map, mu_n = 500) and compare with the
   simulated sign.
2. Prior mass P(r > r*) for tau_p log-uniform on [1, 100] ns and several priors on r, for mu_n = 500 and 50 (crossover from the mu_n
   design map).
Inputs : results/paper_design_map_merged/rstar.csv, results/paper_design_map_mu/rstar_mu.csv, results/paper_uncertainty_ns_merged/ns_samples_all.csv
Output : results/paper_prior_sensitivity/prior_sensitivity.csv, prior_check.csv (printed summary)
"""
import os, csv, math
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = lambda *p: os.path.join(ROOT, "results", *p)
OUT = R("paper_prior_sensitivity"); os.makedirs(OUT, exist_ok=True)
CONDS = ["AM1.5G", "6500K_1000lx"]


def load_rstar():
    d = {}
    for row in csv.DictReader(open(R("paper_design_map_merged", "rstar.csv"))):
        if float(row["Sback"]) != 1e5:
            continue
        d.setdefault((500, row["Condition"]), []).append((float(row["tauP_ns"]), float(row["r_star"])))
    for row in csv.DictReader(open(R("paper_design_map_mu", "rstar_mu.csv"))):
        d.setdefault((int(float(row["muN"])), row["Condition"]), []).append((float(row["tauP_ns"]), float(row["r_star"])))
    return {k: np.array(sorted(v)) for k, v in d.items()}


def rstar_fun(tab):
    lt, lr = np.log(tab[:, 0]), np.log(tab[:, 1])
    return lambda tau_ns: np.exp(np.interp(np.log(tau_ns), lt, lr))   # clamps outside the grid


rs = load_rstar()

# ---- 1. check against the ns sampling ----
rows = list(csv.DictReader(open(R("paper_uncertainty_ns_merged", "ns_samples_all.csv"))))
pcs = {}
for r in rows:
    pcs.setdefault((r["Sample"], r["Condition"]), {})[r["Architecture"]] = r
check = []
for cond in CONDS:
    f = rstar_fun(rs[(500, cond)])
    n = agree = pred_rhj = obs_rhj = 0; margin_dis = []
    for (s, c), a in pcs.items():
        if c != cond or "RBSF-GaP" not in a or "RHJ" not in a:
            continue
        p_rb, p_rhj = float(a["RBSF-GaP"]["PCE_pct"]), float(a["RHJ"]["PCE_pct"])
        if not (math.isfinite(p_rb) and math.isfinite(p_rhj)):
            continue
        tau_ns = float(a["RHJ"]["tauP"]) * 1e9; ratio = float(a["RHJ"]["ratio"])
        pred = ratio > f(tau_ns); obs = p_rhj > p_rb
        n += 1; agree += pred == obs; pred_rhj += pred; obs_rhj += obs
        if pred != obs:
            margin_dis.append(abs(p_rhj - p_rb))
    check.append(dict(Condition=cond, n=n, agree=agree, predicted_RHJ_ahead=pred_rhj, observed_RHJ_ahead=obs_rhj,
                      max_margin_disagreeing=max(margin_dis) if margin_dis else 0.0))
with open(os.path.join(OUT, "prior_check.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(check[0])); w.writeheader(); w.writerows(check)

# ---- 2. prior mass ----
rng = np.random.default_rng(20261107)
M = 2_000_000
tau = np.exp(rng.uniform(math.log(1), math.log(100), M))          # tau_p log-uniform 1-100 ns (as the ns sampling)


def logu(lo, hi):
    return np.exp(rng.uniform(math.log(lo), math.log(hi), M))


def logtri(lo, mode, hi):
    return np.exp(rng.triangular(math.log(lo), math.log(mode), math.log(hi), M))


priors = {
    "log-uniform [1, 30] (sampling)": logu(1, 30),
    "log-uniform [1, 10]": logu(1, 10),
    "log-uniform [1, 100]": logu(1, 100),
    "log-uniform [1, 300]": logu(1, 300),
    "uniform [1, 30]": rng.uniform(1, 30, M),
    "log-triangular [1, 1, 100] (mode 1)": logtri(1, 1, 100),
    "log-triangular [1, 10, 100] (mode 10)": logtri(1, 10, 100),
    "log-triangular [1, 100, 100] (mode 100)": logtri(1, 100, 100),
}
res = []
for mu in (500, 150, 100, 50):
    for cond in CONDS:
        if (mu, cond) not in rs:
            continue
        rstar = rstar_fun(rs[(mu, cond)])(tau)
        for name, r in priors.items():
            p = float(np.mean(r > rstar))
            res.append(dict(muN=mu, Condition=cond, prior_r=name, P_RHJ_ahead=round(p, 4), P_RBSF_ahead=round(1 - p, 4)))
with open(os.path.join(OUT, "prior_sensitivity.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(res[0])); w.writeheader(); w.writerows(res)

for c in check:
    print("CHECK", c)
for x in res:
    print("%4d %-13s %-42s RHJ ahead %.3f" % (x["muN"], x["Condition"], x["prior_r"], x["P_RHJ_ahead"]))
