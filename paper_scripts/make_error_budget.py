"""Error budget of the key PCE differences between architectures (round-5 review R1-8), nominal devices (0.6 us), AM1.5G and 6500 K/1000 lx.
Each entry is the change of the DIFFERENCE (points) caused by one numerical or modelling choice:
  mesh 4x       : D(4x mesh) - D(nominal)                       results/paper_mesh_richardson/richardson_summary.csv
  mesh extrap.  : D(Richardson extrapolation) - D(nominal), where the order is defined for both curves
  interface mesh: largest |D(refined heterointerfaces) - D(nominal)| results/paper_interface_mesh/ (chi_GaP = 3.8)
  voltage step  : D(0.01 V) - D(0.02 V), unified extractor        results/paper_mesh_richardson/unified_step_check.csv (1000 lx only)
  cross-solver  : D(Sesame) - D(MATLAB), SRH-only devices         Supplementary Table tab:supp-sesame (values below)
  Fermi-Dirac   : D(F_1/2) - D(Boltzmann)                         results/paper_fd12_jv/fd12_jv_summary.csv
  Bernoulli     : bound 2 x 2e-3 from the expm1 recomputation      results/paper_r4_analyses/bernoulli_fix_comparison.csv
Output: results/paper_r4_analyses/error_budget.csv and latex/tables/supp_error_budget.tex"""
import glob, os
import numpy as np, pandas as pd

root = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); R = os.path.join(root, "results")
PAIRS = [("RBSF-GaP", "RHJ"), ("RHJ", "F-HJ"), ("RBSF-GaP", "F-HJ"), ("RBSF-GaP", "HJ-0")]
CONDS = ["AM1.5G", "6500K_1000lx"]
# Sesame cross-check, Delta PCE = Sesame - MATLAB (points), Supplementary Table tab:supp-sesame (nominal, SRH only)
SES = {"AM1.5G": {"HJ-0": 0.045, "F-HJ": 0.049, "RBSF-GaP": 0.050, "RHJ": 0.050},
       "6500K_1000lx": {"F-HJ": 0.024, "RBSF-GaP": 0.024, "RHJ": 0.022}}

ms = pd.read_csv(os.path.join(R, "paper_mesh_richardson", "richardson_summary.csv")); ms = ms[ms.metric == "PCE"].set_index(["Architecture", "Condition"])
im = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(os.path.join(R, "paper_interface_mesh", "interface_mesh_*of6.csv")))], ignore_index=True)
im = im[(np.isclose(im.tauP_s, 6e-7)) & (np.isclose(im.chiGaP, 3.8))]
st = pd.read_csv(os.path.join(R, "paper_mesh_richardson", "unified_step_check.csv")); st = st[st.Spectrum == "6500K"].set_index("Architecture")
fd = pd.read_csv(os.path.join(R, "paper_fd12_jv", "fd12_jv_summary.csv")).set_index(["Architecture", "Condition"])

rows = []
for cond in CONDS:
    for a, b in PAIRS:
        D0 = ms.loc[(a, cond), "nominal"] - ms.loc[(b, cond), "nominal"]
        dm4 = (ms.loc[(a, cond), "x4"] - ms.loc[(b, cond), "x4"]) - D0
        ea, eb = ms.loc[(a, cond), "extrapolated"], ms.loc[(b, cond), "extrapolated"]
        dme = (ea - eb) - D0 if np.isfinite(ea) and np.isfinite(eb) else np.nan
        g = im[im.Condition == cond]
        def dI(s):
            pa = g[(g.Architecture == a) & (g.SlabNodes == s)].PCE_pct; pb = g[(g.Architecture == b) & (g.SlabNodes == s)].PCE_pct
            return float(pa.iloc[0] - pb.iloc[0]) if len(pa) and len(pb) else np.nan
        base = dI(0); ref = [dI(s) - base for s in sorted(g.SlabNodes.unique()) if s > 0]
        dif = max(ref, key=abs) if ref and np.all(np.isfinite(ref)) else np.nan
        dst = (st.loc[a, "unified_0.01V"] - st.loc[b, "unified_0.01V"]) - (st.loc[a, "unified_0.02V"] - st.loc[b, "unified_0.02V"]) if cond != "AM1.5G" else np.nan
        dse = SES[cond][a] - SES[cond][b] if a in SES[cond] and b in SES[cond] else np.nan
        dfd = (fd.loc[(a, cond), "PCEFD"] - fd.loc[(b, cond), "PCEFD"]) - (fd.loc[(a, cond), "PCEB"] - fd.loc[(b, cond), "PCEB"])
        terms = [dm4, dif, dst, dse, dfd, 0.004]
        rss = float(np.sqrt(np.nansum(np.square(terms))))
        rows.append((cond, "%s - %s" % (a, b), D0, dm4, dme, dif, dst, dse, dfd, 0.004, rss))
T = pd.DataFrame(rows, columns=["Condition", "Difference", "D_nominal", "mesh4x", "mesh_extrap", "interface_mesh", "voltage_step", "cross_solver", "fermi_dirac", "bernoulli_bound", "rss"])
T.to_csv(os.path.join(R, "paper_r4_analyses", "error_budget.csv"), index=False, lineterminator=chr(10), float_format="%.4f")
pd.set_option("display.width", 220); print(T.round(3).to_string(index=False))

f = lambda v: "--" if not np.isfinite(v) else ("%+.3f" % v).replace("-", "$-$")
nm = {"RBSF-GaP": "RBSF--GaP", "RHJ": "RHJ", "F-HJ": "F-HJ", "HJ-0": "HJ-0"}
L = [r"\begin{table}[htbp]", r"\centering",
     r"\caption{Error budget of the key PCE differences (points) at the nominal parameters (0.6~$\mu$s). Each column gives the change of the difference caused by one choice: four-times denser mesh; refined heterointerfaces (largest change); voltage step 0.01 instead of 0.02~V (1000~lx only); Sesame instead of the present solver (SRH-only devices); $F_{1/2}$ Fermi--Dirac instead of Boltzmann statistics; and the bound from the rounding error of the Bernoulli function. RSS: root sum of squares. In the nanosecond regime the mesh term reaches 0.10 points (Section~\ref{sec:supp-r4B}).}",
     r"\label{tab:supp-errbudget}", r"\footnotesize", r"\setlength{\tabcolsep}{3.5pt}",
     r"\begin{tabular}{llrrrrrrrr}", r"\hline",
     r"Condition & Difference & Nominal & Mesh $4\times$ & Interface & $\Delta V$ & Sesame & $F_{1/2}$ & Bernoulli & RSS \\", r"\hline"]
for i, r in T.iterrows():
    a, b = r.Difference.split(" - ")
    c = ("AM1.5G" if r.Condition == "AM1.5G" else r"1000~lx") if i % len(PAIRS) == 0 else ""
    if i == len(PAIRS): L.append(r"\hline")
    L.append("%s & %s$-$%s & %s & %s & %s & %s & %s & %s & $\\le$0.004 & %.3f \\\\" % (c, nm[a], nm[b], f(r.D_nominal), f(r.mesh4x), f(r.interface_mesh), f(r.voltage_step), f(r.cross_solver), f(r.fermi_dirac), r.rss))
L += [r"\hline", r"\end{tabular}", r"\end{table}"]
open(os.path.join(root, "latex", "tables", "supp_error_budget.tex"), "w", encoding="utf-8", newline=chr(10)).write(chr(10).join(L) + chr(10))
