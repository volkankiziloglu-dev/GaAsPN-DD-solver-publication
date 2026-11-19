"""Analysis of the round-4 group-B runs (docs/HAKEM_SIMULASYONU_2026-11-07.md, roadmap B6-B11).
Inputs: results/paper_r4_analyses/{L200,OPT,SRH}_*of*.csv (cases/run_paper_r4_analyses.m), results/paper_mesh_richardson/<code>_tau<t>_x<s>_<cond>.csv
(cases/run_paper_mesh_richardson.m with tauNs = 1, 10), results/paper_offset_thermionic/offset_thermionic_c4{22,3}_*.csv, results/paper_fd_check5/
(cases/run_paper_fd_check5.m) with results/paper_fd_check4/ for the Boltzmann reference, results/paper_design_map_merged/rstar.csv.
Outputs (results/paper_r4_analyses/): rstar_200lx.csv, sback_200lx.csv, opt_summary.csv, srh_summary.csv, mesh_tau_summary.csv, te_highchi.csv,
fd5_summary.csv; everything is also printed. Sections are skipped when their inputs are missing."""
import glob, os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from unified_metrics import pv_metrics

root = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); R = os.path.join(root, "results"); OUT = os.path.join(R, "paper_r4_analyses")
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40)


def load(part):
    fs = sorted(glob.glob(os.path.join(OUT, f"{part}_*of*.csv")))
    return pd.concat([pd.read_csv(f) for f in fs], ignore_index=True) if fs else None


def rstar_from(rb_pce, rh):   # rh: DataFrame with ratio, PCE_pct (RHJ); linear in log10 r, as analyze_design_map.py
    rh = rh.sort_values("ratio"); r = rh.ratio.values; d = rh.PCE_pct.values - rb_pce; ok = np.isfinite(d); r, d = r[ok], d[ok]
    if len(r) < 2: return np.nan
    if d[0] >= 0: return 1.0
    if d[-1] < 0: return np.inf
    k = np.where(d >= 0)[0][0]; lr = np.log10(r)
    return 10 ** (lr[k - 1] + (0 - d[k - 1]) * (lr[k] - lr[k - 1]) / (d[k] - d[k - 1]))


RS = pd.read_csv(os.path.join(R, "paper_design_map_merged", "rstar.csv")); RS = RS[RS.Sback == 1e5]

# ---------- B7: 200 lx ----------
L = load("L200")
if L is not None:
    L = L[L.Condition == "6500K_200lx"]
    rows = []
    for tp, g in L[L.Sback == 1e5].groupby("tauP_ns"):
        rb = g[(g.Architecture == "RBSF-GaP")]
        if rb.empty: continue
        rs = rstar_from(float(rb.PCE_pct.iloc[0]), g[g.Architecture == "RHJ"])
        r1000 = RS[(RS.Condition == "6500K_1000lx") & (RS.tauP_ns == tp)].r_star; rAM = RS[(RS.Condition == "AM1.5G") & (RS.tauP_ns == tp)].r_star
        rh1 = g[(g.Architecture == "RHJ") & (g.ratio == 1)].PCE_pct
        rows.append((tp, float(rb.PCE_pct.iloc[0]), float(rh1.iloc[0]) if len(rh1) else np.nan, rs, float(r1000.iloc[0]) if len(r1000) else np.nan, float(rAM.iloc[0]) if len(rAM) else np.nan))
    A = pd.DataFrame(rows, columns=["tauP_ns", "PCE_RBSF_200lx", "PCE_RHJ_r1_200lx", "r_star_200lx", "r_star_1000lx", "r_star_AM"])
    A.to_csv(os.path.join(OUT, "rstar_200lx.csv"), index=False, float_format="%.3f"); print("== B7 r* at 200 lx"); print(A.round(3).to_string(index=False))
    sb = L[(L.ratio == 1) & L.Architecture.isin(["HJ-0", "HJ-U", "F-HJ", "RBSF-GaP"])].pivot_table(index=["tauP_ns", "Sback"], columns="Architecture", values="PCE_pct")
    for a in ["HJ-0", "HJ-U", "F-HJ"]: sb["RBSF-" + a] = sb["RBSF-GaP"] - sb[a]
    sb.to_csv(os.path.join(OUT, "sback_200lx.csv"), float_format="%.3f"); print("== B7 rear-contact controls at 200 lx"); print(sb.round(2).to_string())
    print("classes:", L.Class.value_counts().to_dict())

# ---------- B8: thickness / doping ----------
O = load("OPT"); O2 = load("OPT2")
if O is not None and O2 is not None: O = pd.concat([O, O2], ignore_index=True)
if O is not None:
    rows = []
    for (tp, cond, a), g in O.groupby(["tauP_ns", "Condition", "Architecture"]):
        g = g[np.isfinite(g.PCE_pct)]
        if g.empty: continue
        b = g.loc[g.PCE_pct.idxmax()]; nomr = g[(g.tScale == 1) & (np.isclose(g.Nabs, 1e17, rtol=1e-6, atol=0))]
        rows.append((tp, cond, a, float(nomr.PCE_pct.iloc[0]) if len(nomr) else np.nan, b.PCE_pct, b.tScale * 1.5, b.Nabs, len(g)))
    S = pd.DataFrame(rows, columns=["tauP_ns", "Condition", "Architecture", "PCE_nominal_stack", "PCE_best", "best_thickness_um", "best_Nabs", "n"])
    S["gain"] = S.PCE_best - S.PCE_nominal_stack
    S.to_csv(os.path.join(OUT, "opt_summary.csv"), index=False, float_format="%.3f"); print("== B8 per-architecture optimum (grid)"); print(S.round(3).to_string(index=False))
    for (tp, cond), g in S.groupby(["tauP_ns", "Condition"]):
        print(tp, cond, "nominal order:", " > ".join(g.sort_values("PCE_nominal_stack", ascending=False).Architecture), "| optimized order:", " > ".join(g.sort_values("PCE_best", ascending=False).Architecture))
    print("errors/NaN:", int((~np.isfinite(O.PCE_pct)).sum()))

# ---------- B11: injection-dependent SRH ----------
Sr = load("SRH")
if Sr is not None:
    P = Sr.pivot_table(index=["tauP_ns", "ratio", "EtRelEi_eV", "Condition"], columns="Architecture", values="PCE_pct").reset_index()
    P["RBSF_minus_RHJ"] = P["RBSF-GaP"] - P["RHJ"]
    L2 = pd.read_csv(os.path.join(OUT, "rstar_200lx.csv")) if os.path.exists(os.path.join(OUT, "rstar_200lx.csv")) else None
    def rmap(cond, tau):
        if cond == "6500K_200lx":
            if L2 is None: return np.nan
            t = L2.sort_values("tauP_ns"); return float(np.exp(np.interp(np.log(tau), np.log(t.tauP_ns), np.log(np.minimum(t.r_star_200lx, 1e4)))))
        t = RS[RS.Condition == cond].sort_values("tauP_ns"); return float(np.exp(np.interp(np.log(tau), np.log(t.tauP_ns), np.log(np.minimum(t.r_star, 1e4)))))
    P["r_star_map"] = [rmap(c, t) for c, t in zip(P.Condition, P.tauP_ns)]
    P["map_predicts_RHJ"] = P.ratio > P.r_star_map; P["RHJ_ahead"] = P.RBSF_minus_RHJ < 0
    P.to_csv(os.path.join(OUT, "srh_summary.csv"), index=False, float_format="%.3f"); print("== B11 single SRH centre in both absorbers"); print(P.round(3).to_string(index=False))
    print("sign predicted by constant-r map:", int((P.map_predicts_RHJ == P.RHJ_ahead).sum()), "of", len(P))

# ---------- B11b: single-centre ratio scan ----------
S2 = load("SRH2")
if S2 is not None:
    P2 = S2.pivot_table(index=["tauP_ns", "ratio", "EtRelEi_eV", "Condition"], columns="Architecture", values="PCE_pct").reset_index()
    P2["RHJ_minus_RBSF"] = P2["RHJ"] - P2["RBSF-GaP"]
    P2.to_csv(os.path.join(OUT, "srh2_summary.csv"), index=False, float_format="%.3f"); print("== B11b single-centre ratio scan"); print(P2.round(3).to_string(index=False))
    rows = []
    for (tn, cond), g in P2[P2.EtRelEi_eV == 0].groupby(["tauP_ns", "Condition"]):
        g = g.sort_values("ratio"); r = g.ratio.values; d = g.RHJ_minus_RBSF.values
        if d[0] >= 0: rs = 1.0
        elif d[-1] < 0: rs = np.inf
        else:
            k = np.where(d >= 0)[0][0]; lr = np.log10(r); rs = 10 ** (lr[k - 1] - d[k - 1] * (lr[k] - lr[k - 1]) / (d[k] - d[k - 1]))
        cmap = "6500K_1000lx" if cond == "6500K_1000lx" else cond
        if cond == "6500K_200lx":
            t = pd.read_csv(os.path.join(OUT, "rstar_200lx.csv")); rm = float(np.exp(np.interp(np.log(tn), np.log(t.tauP_ns), np.log(t.r_star_200lx))))
        else:
            t = RS[RS.Condition == cond].sort_values("tauP_ns"); rm = float(np.exp(np.interp(np.log(tn), np.log(t.tauP_ns), np.log(t.r_star))))
        rows.append((tn, cond, rs, rm))
    X = pd.DataFrame(rows, columns=["tau_n0_ns", "Condition", "r_star_single_centre", "r_star_map"])
    X.to_csv(os.path.join(OUT, "srh2_rstar.csv"), index=False, float_format="%.3f"); print(X.round(2).to_string(index=False))

# ---------- B6a: mesh at 1 and 10 ns ----------
d = os.path.join(R, "paper_mesh_richardson")
Sm = pd.read_csv(os.path.join(R, "paper_indoor_production", "indoor_production_summary.csv"))
Pin = {"6500K_1000lx": float(Sm[(Sm.Architecture == "RBSF-GaP") & (Sm.Spectrum == "6500K") & (Sm.Lux == 1000)].Pin_Wm2.iloc[0]), "AM1.5G": 1000.3706555734421}
rows = []
for tau in (1, 10):
    for a in ["HJ-0", "F-HJ", "RBSF-GaP", "RHJ"]:
        for cond in Pin:
            for sc in (1, 2, 4):
                f = os.path.join(d, f"{a}_tau{tau}_x{sc}_{cond}.csv")
                if not os.path.exists(f): continue
                c = pd.read_csv(f); c = c[np.isfinite(c.J_Acm2)]; m = pv_metrics(c.V_V, c.J_Acm2, Pin[cond])
                rows.append((tau, a, cond, sc, m["PCE"], m["Voc"] * 1e3, int(c.Accepted.sum()), len(c)))
if rows:
    M = pd.DataFrame(rows, columns=["tau_ns", "Architecture", "Condition", "scale", "PCE", "Voc_mV", "nAccepted", "nPoints"])
    M.to_csv(os.path.join(OUT, "mesh_tau_points.csv"), index=False, float_format="%.4f")
    W = M.pivot_table(index=["tau_ns", "Condition", "scale"], columns="Architecture", values="PCE")
    for a, b in [("RBSF-GaP", "RHJ"), ("RHJ", "F-HJ"), ("RBSF-GaP", "F-HJ"), ("RBSF-GaP", "HJ-0")]:
        if a in W and b in W: W[a + " - " + b] = W[a] - W[b]
    W.to_csv(os.path.join(OUT, "mesh_tau_summary.csv"), float_format="%.4f"); print("== B6 mesh at 1 and 10 ns (PCE)"); print(W.round(3).to_string())

# ---------- B6b: thermionic emission at high chi_GaP ----------
fs = sorted(glob.glob(os.path.join(R, "paper_offset_thermionic", "offset_thermionic_c4*_*.csv")))
if fs:
    T = pd.concat([pd.read_csv(f) for f in fs]); T.to_csv(os.path.join(OUT, "te_highchi.csv"), index=False)
    T["model"] = np.where(T.v_cm_s.isna(), "SG", "TE " + T.v_cm_s.astype(str))
    print("== B6 thermionic emission at high chi_GaP"); print(T.pivot_table(index=["chiGaP", "Condition", "model"], columns="Architecture", values="PCE_pct").round(3).to_string())

# ---------- B10: consistent degenerate SG ----------
f5 = sorted(glob.glob(os.path.join(R, "paper_fd_check5", "fd5_*_[12].csv")))
if f5:
    F = pd.concat([pd.read_csv(f) for f in f5]); F4 = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(os.path.join(R, "paper_fd_check4", "fd4_*_[12].csv")))])
    F = pd.concat([F, F4[F4.Variant != "Boltzmann"]])
    W = F.pivot_table(index=["Condition", "Variant"], columns="Architecture", values="PCE_pct")
    for a, b in [("RBSF-GaP", "RHJ"), ("RHJ", "F-HJ"), ("RBSF-GaP", "HJ-0")]: W[a + " - " + b] = W[a] - W[b]
    W.to_csv(os.path.join(OUT, "fd5_summary.csv"), float_format="%.4f"); print("== B10 Fermi-Dirac discretizations (PCE)"); print(W.round(3).to_string())
    print(F[["Architecture", "Condition", "Variant", "Voc_mV", "Status"]].to_string(index=False))
