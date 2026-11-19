"""Dependence of the crossover ratio r*(tau_p) on the GaP electron affinity chi_GaP (simulated referee round 4, W3/R2/R4,
docs/HAKEM_SIMULASYONU_2026-11-07.md). chi_abs = 4.0 eV, so Delta E_c(GaP/absorber) = chi_abs - chi_GaP: 3.7 -> +0.30 eV (type I),
3.8 (nominal) -> +0.20, 4.0 -> 0, 4.15 -> -0.15, 4.3 -> -0.30 eV (type II, GaP conduction band below the absorber's).
Inputs : results/paper_design_map_chi/design_map_chi<chi>_{a,b}.csv (RBSF-GaP and RHJ, S_back = 1e5, mu_n = 500),
         results/paper_design_map_merged/design_map_[abcd].csv (nominal chi_GaP = 3.8),
         results/paper_uncertainty_ns_merged/ns_samples_all.csv (nanosecond sampling, 60 LHS samples).
1. r*(tau_p, chi_GaP) table (linear interpolation of PCE_RHJ - PCE_RBSF versus log10 r, as analyze_design_map.py; r* = 1 is reported
   when RHJ is already ahead at r = 1 (bound), inf when RHJ never catches up within r <= 60).
2. Re-classification of the ns sampling with a chi-aware map (log r* interpolated linearly in chi_GaP and log tau_p) versus the
   chi-blind nominal map.
3. Prior mass P(RHJ ahead) for tau_p log-uniform [1, 100] ns, r log-uniform [1, 30] and chi_GaP uniform on [3.7, 4.3] (sampling prior),
   [3.7, 4.0] (type-I only) or fixed at 3.8; seed 20261108, M = 2e6.
Output : results/paper_design_map_chi/rstar_chi.csv, chi_check.csv, chi_prior.csv (printed summary)."""
import glob, os, math
import numpy as np, pandas as pd

root = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); R = os.path.join(root, "results")
OUT = os.path.join(R, "paper_design_map_chi")
CONDS = ["AM1.5G", "6500K_1000lx"]
CAP = 1e4   # stand-in for "no crossover within r <= 60" in log interpolation
frames = []
for f in sorted(glob.glob(os.path.join(OUT, "design_map_chi*_[ab].csv"))):
    chi = float(os.path.basename(f)[len("design_map_chi"):-6]); d = pd.read_csv(f); d["chiGaP"] = chi; frames.append(d)
m = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(os.path.join(R, "paper_design_map_merged", "design_map_[abcd].csv")))])
m = m[(m.Sback == 1e5) & m.Architecture.isin(["RBSF-GaP", "RHJ"])].copy(); m["chiGaP"] = 3.8; frames.append(m)
D = pd.concat(frames, ignore_index=True)

rows = []
for (chi, cond, tp), g in D.groupby(["chiGaP", "Condition", "tauP_ns"]):
    rb = g[g.Architecture == "RBSF-GaP"]
    if rb.empty: continue
    rb = rb.iloc[0]; rh = g[g.Architecture == "RHJ"].sort_values("ratio")
    r = rh.ratio.values; diff = (rh.PCE_pct - rb.PCE_pct).values; ok = np.isfinite(diff); r, diff = r[ok], diff[ok]
    if len(r) < 2: continue
    lr = np.log10(r)
    if diff[0] >= 0: rstar = 1.0
    elif diff[-1] < 0: rstar = np.inf
    else:
        k = np.where(diff >= 0)[0][0]; rstar = 10 ** (lr[k - 1] + (0 - diff[k - 1]) * (lr[k] - lr[k - 1]) / (diff[k] - diff[k - 1]))
    rows.append((chi, cond, tp, rb.PCE_pct, diff[0], rstar, int((~ok).sum())))
RS = pd.DataFrame(rows, columns=["chiGaP", "Condition", "tauP_ns", "PCE_RBSF", "RHJ_minus_RBSF_r1", "r_star", "n_nonfinite"])
RS.to_csv(os.path.join(OUT, "rstar_chi.csv"), index=False, float_format="%.3f")
pd.set_option("display.width", 200)
for cond in CONDS:
    s = RS[RS.Condition == cond]
    print(cond, "r*"); print(s.pivot(index="tauP_ns", columns="chiGaP", values="r_star").round(2).to_string())
    print(cond, "RHJ - RBSF at r = 1 (points)"); print(s.pivot(index="tauP_ns", columns="chiGaP", values="RHJ_minus_RBSF_r1").round(2).to_string())
    print(cond, "PCE RBSF-GaP"); print(s.pivot(index="tauP_ns", columns="chiGaP", values="PCE_RBSF").round(2).to_string())


def rstar_fn(cond, chi_aware):
    s = RS[RS.Condition == cond]
    P = s.pivot(index="tauP_ns", columns="chiGaP", values="r_star").sort_index()
    lt = np.log(P.index.values.astype(float)); chis = P.columns.values.astype(float)
    L = np.log(np.minimum(P.values.astype(float), CAP))
    if not chi_aware:
        j = list(chis).index(3.8)
        return lambda tau_ns, chi: np.exp(np.interp(np.log(tau_ns), lt, L[:, j]))
    def f(tau_ns, chi):
        tau_ns = np.atleast_1d(tau_ns); chi = np.atleast_1d(chi) * np.ones_like(tau_ns)
        Lt = np.array([np.interp(np.log(tau_ns), lt, L[:, j]) for j in range(len(chis))])   # (nchi, n)
        out = np.empty(tau_ns.shape)
        for i in range(len(tau_ns)):
            out[i] = np.interp(chi[i], chis, Lt[:, i])
        return np.exp(out)
    return f


S = pd.read_csv(os.path.join(R, "paper_uncertainty_ns_merged", "ns_samples_all.csv"))
chk = []
for cond in CONDS:
    t = S[S.Condition == cond]
    P = t.pivot(index="Sample", columns="Architecture", values="PCE_pct")
    prm = t.drop_duplicates("Sample").set_index("Sample").loc[P.index]
    ok = np.isfinite(P["RBSF-GaP"]) & np.isfinite(P["RHJ"])
    P, prm = P[ok], prm[ok]
    obs = (P["RHJ"] > P["RBSF-GaP"]).values; margin = (P["RHJ"] - P["RBSF-GaP"]).values
    tau = prm.tauP.values * 1e9; r = prm.ratio.values; chi = prm.chiGaP.values
    for aware in (False, True):
        pred = r > rstar_fn(cond, aware)(tau, chi)
        bad = pred != obs
        chk.append(dict(Condition=cond, map="chi-aware" if aware else "nominal chi 3.8", n=int(len(obs)), agree=int((~bad).sum()),
                        predicted_RHJ_ahead=int(pred.sum()), observed_RHJ_ahead=int(obs.sum()),
                        max_margin_disagreeing=float(np.abs(margin[bad]).max()) if bad.any() else 0.0,
                        chi_of_disagreeing=" ".join("%.2f" % x for x in chi[bad]), samples_disagreeing=" ".join(str(x) for x in P.index[bad])))
    for lo, hi, lab in ((3.7, 4.3, "all"), (3.7, 4.0, "type-I chi<=4.0"), (4.0, 4.3, "type-II chi>4.0")):
        sel = (chi >= lo) & (chi <= hi)
        chk.append(dict(Condition=cond, map="observed, " + lab, n=int(sel.sum()), agree=np.nan, predicted_RHJ_ahead=np.nan,
                        observed_RHJ_ahead=int(obs[sel].sum()), max_margin_disagreeing=np.nan, chi_of_disagreeing="", samples_disagreeing=""))
C = pd.DataFrame(chk); C.to_csv(os.path.join(OUT, "chi_check.csv"), index=False, float_format="%.3f")
print(C.to_string())

rng = np.random.default_rng(20261108); M = 2_000_000
tau = np.exp(rng.uniform(math.log(1), math.log(100), M)); r = np.exp(rng.uniform(math.log(1), math.log(30), M)); u = rng.uniform(0, 1, M)
res = []
for cond in CONDS:
    fa = rstar_fn(cond, True)
    # vectorised chi-aware evaluation
    s = RS[RS.Condition == cond]; Pm = s.pivot(index="tauP_ns", columns="chiGaP", values="r_star").sort_index()
    lt = np.log(Pm.index.values.astype(float)); chis = Pm.columns.values.astype(float); L = np.log(np.minimum(Pm.values.astype(float), CAP))
    Lt = np.array([np.interp(np.log(tau), lt, L[:, j]) for j in range(len(chis))])
    for lo, hi, lab in ((3.7, 4.3, "uniform [3.7, 4.3] (sampling prior)"), (3.7, 4.0, "uniform [3.7, 4.0] (type I)"), (3.8, 3.8, "fixed 3.8 (nominal)"),
                        (4.0, 4.3, "uniform [4.0, 4.3] (type II)")):
        chi = lo + (hi - lo) * u
        j = np.clip(np.searchsorted(chis, chi) - 1, 0, len(chis) - 2); w = (chi - chis[j]) / (chis[j + 1] - chis[j])
        lr = (1 - w) * Lt[j, np.arange(M)] + w * Lt[j + 1, np.arange(M)]
        res.append(dict(Condition=cond, chi_prior=lab, P_RHJ_ahead=round(float(np.mean(np.log(r) > lr)), 4)))
Q = pd.DataFrame(res); Q.to_csv(os.path.join(OUT, "chi_prior.csv"), index=False)
print(Q.to_string())
