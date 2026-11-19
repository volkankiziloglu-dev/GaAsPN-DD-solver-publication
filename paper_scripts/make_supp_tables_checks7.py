"""LaTeX tables of the checks requested by the third simulated review (docs/HAKEM_SIMULASYONU_2026-11-06.md), read from the result CSVs.
Writes latex/tables/supp_checks7_<name>.tex (one file per table), included by latex/supp_mat_en.tex.
Inputs: results/paper_sfront, results/paper_pin (parts U, V), results/paper_referee_checks_merged/referee_checks_A.csv, results/paper_pin (N, M),
        results/paper_prior_sensitivity, external_validation/scaps_output/scaps_prod_crosscheck.csv, results/paper_mesh_cvgen,
        results/paper_interface_mesh, results/paper_fd_check4 (if complete), results/paper_design_map_chi (round-4 review, 2026-11-07).
"""
import os, csv, glob, math

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = lambda *p: os.path.join(ROOT, *p)
OUT = R("latex", "tables"); os.makedirs(OUT, exist_ok=True)
ARCH = ["HJ-0", "F-HJ", "RBSF-GaP", "RHJ"]
TEX = {"HJ-0": "HJ-0", "F-HJ": "F-HJ", "RBSF-GaP": "RBSF--GaP", "RHJ": "RHJ", "HJ-U": "HJ-U", "PIN": "p--i--n"}


def rows(pattern):
    return [r for f in sorted(glob.glob(R(*pattern))) for r in csv.DictReader(open(f))]


def sgn(x):
    v = float(x)
    return ("$-$%.2f" % abs(v)) if v < 0 and round(v, 2) != 0 else "%.2f" % abs(v)


def write(name, lines):
    open(os.path.join(OUT, "supp_checks7_%s.tex" % name), "w", encoding="utf-8", newline=chr(10)).write("\n".join(lines) + "\n")


# ---------- 1. front-contact velocity ----------
sf = {(r["Architecture"], float(r["tauP_s"]), r["Sfront"], r["Condition"]): float(r["PCE_pct"]) for r in rows(["results", "paper_sfront", "sfront_*of6.csv"])}
L = ["\\begin{table*}[htbp]", "\\centering",
     "\\caption{PCE (\\%) versus the front-contact recombination velocity $S_{\\rm front}$ of the minority holes in the n$^+$-GaP window ($\\infty$: ideal ohmic contact of the main calculations); AM1.5G / 6500~K, 1000~lx; nominal parameters except $\\tau_p$ ($r=1$).}",
     "\\label{tab:supp-sfront}", "\\scriptsize", "\\begin{tabular}{llcccc}", "\\toprule",
     "$\\tau_p$ & $S_{\\rm front}$ (cm\\,s$^{-1}$) & HJ-0 & F-HJ & RBSF--GaP & RHJ \\\\", "\\midrule"]
for tau, tl in [(6e-7, "0.6~$\\mu$s"), (1e-8, "10~ns"), (1e-9, "1~ns")]:
    for S, Sl in [("Inf", "$\\infty$"), ("1000000", "$10^6$"), ("10000", "$10^4$"), ("100", "$10^2$")]:
        cells = ["%.2f / %.2f" % (sf[(a, tau, S, "AM1.5G")], sf[(a, tau, S, "6500K_1000lx")]) for a in ARCH]
        L.append("%s & %s & %s \\\\" % (tl if S == "Inf" else "", Sl, " & ".join(cells)))
    if tau != 1e-9: L.append("\\midrule")
L += ["\\bottomrule", "\\end{tabular}", "\\end{table*}"]
write("sfront", L)

# ---------- 2. HJ-U control ----------
A = {(r["Architecture"], float(r["tauP_s"]), r["Condition"]): float(r["PCE_pct"]) for r in rows(["results", "paper_referee_checks_merged", "referee_checks_A.csv"]) if float(r["ratio"]) == 1}
U = {(float(r["tauP_s"]), r["Condition"], int(float(r["MuN"]))): float(r["PCE_pct"]) for r in rows(["results", "paper_pin", "pin_[UV]_*.csv"])}
M = {(r["Architecture"], float(r["tauP_s"]), r["Condition"], float(r["ratio"])): float(r["PCE_pct"]) for r in rows(["results", "paper_pin", "pin_M_*.csv"])}
L = ["\\begin{table*}[htbp]", "\\centering",
     "\\caption{Homojunction with a p$^+$-GaAsPN rear layer (HJ-U, passivated-rear control) against the four architectures; PCE (\\%) AM1.5G / 6500~K, 1000~lx; $r=1$, $S_{\\rm back}=10^5$~cm\\,s$^{-1}$. Last column: HJ-U minus RBSF--GaP.}",
     "\\label{tab:supp-hju}", "\\scriptsize", "\\begin{tabular}{llcccccc}", "\\toprule",
     "$\\mu_n$ & $\\tau_p$ & HJ-U & HJ-0 & F-HJ & RBSF--GaP & RHJ & HJ-U$-$RBSF--GaP \\\\", "\\midrule"]
for tau, tl in [(6e-7, "0.6~$\\mu$s"), (1e-7, "100~ns"), (1e-8, "10~ns"), (1e-9, "1~ns")]:
    hu = [U[(tau, c, 500)] for c in ("AM1.5G", "6500K_1000lx")]
    oth = [" / ".join("%.2f" % A[(a, tau, c)] for c in ("AM1.5G", "6500K_1000lx")) for a in ARCH]
    d = [hu[i] - A[("RBSF-GaP", tau, c)] for i, c in enumerate(("AM1.5G", "6500K_1000lx"))]
    L.append("%s & %s & %.2f / %.2f & %s & %+.2f / %+.2f \\\\" % ("500" if tau == 6e-7 else "", tl, hu[0], hu[1], " & ".join(oth), d[0], d[1]))
L.append("\\midrule")
for tau, tl in [(6e-7, "0.6~$\\mu$s"), (1e-8, "10~ns"), (1e-9, "1~ns")]:
    hu = [U[(tau, c, 50)] for c in ("AM1.5G", "6500K_1000lx")]
    oth = [" / ".join("%.2f" % M[(a, tau, c, 1.0)] for c in ("AM1.5G", "6500K_1000lx")) for a in ARCH]
    d = [hu[i] - M[("RBSF-GaP", tau, c, 1.0)] for i, c in enumerate(("AM1.5G", "6500K_1000lx"))]
    L.append("%s & %s & %.2f / %.2f & %s & %+.2f / %+.2f \\\\" % ("50" if tau == 6e-7 else "", tl, hu[0], hu[1], " & ".join(oth), d[0], d[1]))
L += ["\\bottomrule", "\\end{tabular}", "\\end{table*}"]
write("hju", L)

# ---------- 3. prior sensitivity ----------
ps = rows(["results", "paper_prior_sensitivity", "prior_sensitivity.csv"])
ck = rows(["results", "paper_prior_sensitivity", "prior_check.csv"])
P = {(int(r["muN"]), r["Condition"], r["prior_r"]): float(r["P_RHJ_ahead"]) for r in ps}
priors = list(dict.fromkeys(r["prior_r"] for r in ps))
L = ["\\begin{table}[htbp]", "\\centering",
     "\\caption{Share of the $(\\tau_p, r)$ prior in which the RHJ leads RBSF--GaP, $P(r>r^*(\\tau_p))$, for $\\tau_p$ log-uniform on 1--100~ns and several priors on $r$; $r^*(\\tau_p)$ from the design map ($S_{\\rm back}=10^5$~cm\\,s$^{-1}$; log--log interpolation; $2\\times10^6$ draws, NumPy PCG64 seed 20261107). AM1.5G / 6500~K, 1000~lx. Check against the nanosecond sampling: the rule $r>r^*$ predicts the simulated sign in %s and %s of 60 samples.}" % (ck[0]["agree"], ck[1]["agree"]),
     "\\label{tab:supp-prior}", "\\scriptsize", "\\begin{tabular}{lcccc}", "\\toprule",
     "Prior on $r$ & $\\mu_n=500$ & $\\mu_n=150$ & $\\mu_n=100$ & $\\mu_n=50$ \\\\", "\\midrule"]
for p in priors:
    name = p.replace("(sampling)", "(sampling)").replace("[", "$[").replace("]", "]$")
    L.append("%s & %s \\\\" % (name, " & ".join("%.2f / %.2f" % (P[(mu, "AM1.5G", p)], P[(mu, "6500K_1000lx", p)]) for mu in (500, 150, 100, 50))))
L += ["\\bottomrule", "\\end{tabular}", "\\end{table}"]
write("prior", L)

# ---------- 4. SCAPS production configuration ----------
sc = rows(["external_validation", "scaps_output", "scaps_prod_crosscheck.csv"])
S = {(r["tag"], r["code"]): r for r in sc}
tags = list(dict.fromkeys(r["tag"] for r in sc))
lab = {"": "AM1.5G, 0.6~$\\mu$s", "_tau10ns": "AM1.5G, 10~ns", "_indoor": "1000~lx", "_indoor200": "200~lx", "_indoor100": "100~lx"}
L = ["\\begin{table*}[htbp]", "\\centering",
     "\\caption{SCAPS-1D versus the present model with the production recombination (midgap SRH, radiative, Auger, interface defects with $S_{\\rm int}=8910$~cm\\,s$^{-1}$ at the GaP/GaAsPN interfaces), identical nodes, parameters and generation profile; one voltage grid (0.02~V) and one extractor. $J_{sc}$ in mA\\,cm$^{-2}$ (indoor: $\\mu$A\\,cm$^{-2}$).}",
     "\\label{tab:supp-scaps-prod}", "\\scriptsize", "\\begin{tabular}{llcccccccc}", "\\toprule",
     " & & \\multicolumn{4}{c}{SCAPS-1D} & \\multicolumn{4}{c}{This model} \\\\",
     "Architecture & Condition & $J_{sc}$ & $V_{oc}$ (mV) & FF (\\%) & PCE (\\%) & $J_{sc}$ & $V_{oc}$ (mV) & FF (\\%) & PCE (\\%) \\\\", "\\midrule"]
for t in tags:
    a = t[2:]
    for k in ("RBSF-GaP", "RHJ", "F-HJ"):
        if a.startswith(k): arch = k; suf = a[len(k):]; break
    f = 1e3 if "indoor" in suf else 1.0
    s, m = S[(t, "SCAPS")], S[(t, "MATLAB")]
    L.append("%s & %s & %.4g & %.1f & %.2f & %.2f & %.4g & %.1f & %.2f & %.2f \\\\" % (
        TEX[arch], lab[suf], float(s["Jsc_mAcm2"]) * f, float(s["Voc_mV"]), float(s["FF_pct"]), float(s["PCE_pct"]),
        float(m["Jsc_mAcm2"]) * f, float(m["Voc_mV"]), float(m["FF_pct"]), float(m["PCE_pct"])))
L += ["\\bottomrule", "\\end{tabular}", "\\end{table*}"]
write("scaps_prod", L)

# ---------- 5. control-volume generation ----------
cv = rows(["results", "paper_mesh_cvgen", "mesh_cvgen_*.csv"])
C = {(r["Architecture"], float(r["MeshScale"]), r["Condition"]): r for r in cv}
rich = {(r["Architecture"], r["Condition"], r["metric"]): r for r in rows(["results", "paper_mesh_richardson", "richardson_summary.csv"])}
L = ["\\begin{table*}[htbp]", "\\centering",
     "\\caption{PCE (\\%) with the generation integrated exactly over each control volume, for mesh scales $1\\times$, $2\\times$, $4\\times$, compared with the nodal generation of the production runs ($1\\times$) and the Richardson extrapolation of the nodal three-level study. AM1.5G / 6500~K, 1000~lx; nominal parameters.}",
     "\\label{tab:supp-cvgen}", "\\scriptsize", "\\begin{tabular}{lccccc}", "\\toprule",
     "Architecture & CV $1\\times$ & CV $2\\times$ & CV $4\\times$ & nodal $1\\times$ & nodal extrapolated \\\\", "\\midrule"]
for a in ARCH:
    cells = []
    for s in (1.0, 2.0, 4.0):
        cells.append("%.3f / %.3f" % (float(C[(a, s, "AM1.5G")]["PCE_pct"]), float(C[(a, s, "6500K_1000lx")]["PCE_pct"])))
    nod = "%.3f / %.3f" % (float(rich[(a, "AM1.5G", "PCE")]["nominal"]), float(rich[(a, "6500K_1000lx", "PCE")]["nominal"]))
    ex = [rich[(a, c, "PCE")]["extrapolated"] for c in ("AM1.5G", "6500K_1000lx")]
    ext = " / ".join(("%.3f" % float(x)) if x not in ("", None) else "--" for x in ex)
    L.append("%s & %s & %s & %s \\\\" % (TEX[a], " & ".join(cells), nod, ext))
L += ["\\bottomrule", "\\end{tabular}", "\\end{table*}"]
write("cvgen", L)

# ---------- 6. interface refinement ----------
im = rows(["results", "paper_interface_mesh", "interface_mesh_*of6.csv"])
I = {(r["Architecture"], float(r["tauP_s"]), float(r["chiGaP"]), int(r["SlabNodes"]), r["Condition"]): float(r["PCE_pct"]) for r in im}
L = ["\\begin{table*}[htbp]", "\\centering",
     "\\caption{PCE (\\%) with local refinement of the heterointerfaces: production mesh / 5~nm slabs with 25 nodes (0.2~nm) / 50 nodes (0.1~nm) on both sides of every GaP/GaAsPN interface. Nominal offset ($\\chi_{\\rm GaP}=3.8$~eV, $\\Delta E_c=\\Delta E_v=0.20$~eV) and conduction-band-dominated offset ($\\chi_{\\rm GaP}=3.65$~eV, 0.35/0.05~eV).}",
     "\\label{tab:supp-ifmesh}", "\\scriptsize", "\\begin{tabular}{lllcccc}", "\\toprule",
     "Condition & $\\tau_p$ & $\\chi_{\\rm GaP}$ & HJ-0 & F-HJ & RBSF--GaP & RHJ \\\\", "\\midrule"]
for c, cl in [("AM1.5G", "AM1.5G"), ("6500K_1000lx", "1000~lx")]:
    for tau, tl in [(6e-7, "0.6~$\\mu$s"), (1e-8, "10~ns")]:
        for chi in (3.8, 3.65):
            cells = ["%.2f / %.2f / %.2f" % tuple(I[(a, tau, chi, n, c)] for n in (0, 25, 50)) for a in ARCH]
            L.append("%s & %s & %.2f & %s \\\\" % (cl, tl, chi, " & ".join(cells)))
L += ["\\bottomrule", "\\end{tabular}", "\\end{table*}"]
write("ifmesh", L)

# ---------- 6b. partial rank correlations and equivalence-margin sensitivity (cases/analyze_sampling_r4_minor.py) ----------
pr = rows(["results", "paper_uncertainty_merged", "r4_prcc.csv"])
if pr:
    PN = {"tauAbs": "$\\tau_{\\rm abs}$", "chiGaP": "$\\chi_{\\rm GaP}$", "Nabs": "$N_{\\rm abs}$", "ratio": "$r$", "tauP": "$\\tau_p$",
          "Sback": "$S_{\\rm back}$", "Sint": "$S_{\\rm int}$"}
    QN = lambda q: q.replace("RBSF-GaP-RHJ", "RBSF--GaP$-$RHJ").replace("PCE RBSF-GaP", "PCE, RBSF--GaP").replace("PCE HJ-0", "PCE, HJ-0").replace("PCE RHJ", "PCE, RHJ")
    keys = list(dict.fromkeys((r["sampling"], r["quantity"], r["parameter"]) for r in pr))
    PR = {(r["sampling"], r["quantity"], r["parameter"], r["condition"]): r for r in pr}
    L = ["\\begin{table}[htbp]", "\\centering",
         "\\caption{Spearman rank correlation and partial rank correlation (PRCC, other sampled inputs partialled out) of the quantities discussed in the main text; AM1.5G / 6500~K, 1000~lx; $n=60$.}",
         "\\label{tab:supp-prcc}", "\\scriptsize", "\\begin{tabular}{lllcc}", "\\toprule",
         "Sampling & Quantity & Input & Spearman & PRCC \\\\", "\\midrule"]
    for s, q, p in keys:
        a, b = PR[(s, q, p, "AM1.5G")], PR[(s, q, p, "6500K_1000lx")]
        L.append("%s & %s & %s & %s / %s & %s / %s \\\\" % (s[-1], QN(q), PN.get(p, p), sgn(a["spearman"]), sgn(b["spearman"]), sgn(a["prcc"]), sgn(b["prcc"])))
    L += ["\\bottomrule", "\\end{tabular}", "\\end{table}"]
    write("prcc", L)
    wm = rows(["results", "paper_uncertainty_merged", "r4_wtl_margins.csv"])
    W = {(r["sampling"], r["condition"], r["A"], r["B"], float(r["margin"])): r for r in wm}
    pairs = list(dict.fromkeys((r["A"], r["B"]) for r in wm))
    L = ["\\begin{table}[htbp]", "\\centering",
         "\\caption{Win/tie/loss counts of A over B for equivalence margins of 0.05, 0.1 (used in the text) and 0.2 points; AM1.5G / 6500~K, 1000~lx.}",
         "\\label{tab:supp-wtl-margin}", "\\scriptsize", "\\setlength{\\tabcolsep}{3pt}", "\\begin{tabular}{llccc}", "\\toprule",
         "Sampling & A $>$ B & 0.05 & 0.1 & 0.2 \\\\", "\\midrule"]
    for s in ("sampling 1", "sampling 2"):
        for a, b in pairs:
            cells = []
            for m in (0.05, 0.1, 0.2):
                x, y = W[(s, "AM1.5G", a, b, m)], W[(s, "6500K_1000lx", a, b, m)]
                cells.append("%s/%s/%s, %s/%s/%s" % (x["wins"], x["ties"], x["losses"], y["wins"], y["ties"], y["losses"]))
            L.append("%s & %s $>$ %s & %s \\\\" % (s[-1], TEX[a], TEX[b], " & ".join(cells)))
        if s == "sampling 1": L.append("\\midrule")
    L += ["\\bottomrule", "\\end{tabular}", "\\end{table}"]
    write("wtl_margin", L)

# ---------- 7. Fermi-Dirac (only if complete) ----------
fd = rows(["results", "paper_fd_check4", "fd4_*_[12].csv"])
if fd and all(math.isfinite(float(r["PCE_pct"])) for r in fd if r["Variant"] == "Boltzmann"):
    F = {(r["Architecture"], int(r["Condition"]), r["Variant"]): r for r in fd}
    L = ["\\begin{table*}[htbp]", "\\centering",
         "\\caption{Boltzmann versus Blakemore (Fermi--Dirac) statistics with full sweeps and the unified extractor; Blakemore with the effective intrinsic density $n_{i,\\rm eff}=\\sqrt{n_{\\rm eq}p_{\\rm eq}}$ in the recombination terms, with the plain SG flux and with the degeneracy-corrected SG flux (generalized Einstein relation). $V_{oc}$ (mV) and PCE (\\%), AM1.5G / 6500~K, 1000~lx.}",
         "\\label{tab:supp-fd4}", "\\scriptsize", "\\begin{tabular}{lccc}", "\\toprule",
         "Architecture & Boltzmann & Blakemore & Blakemore, degenerate SG \\\\", "\\midrule"]
    for a in ARCH:
        cells = []
        for v in ("Boltzmann", "Blakemore", "BlakemoreDegenSG"):
            vals = []
            for c in (1, 2):
                r = F.get((a, c, v))
                vals.append("--" if r is None or not math.isfinite(float(r["PCE_pct"])) else "%.1f, %.2f" % (float(r["Voc_mV"]), float(r["PCE_pct"])))
            cells.append(" / ".join(vals))
        L.append("%s & %s \\\\" % (TEX[a], " & ".join(cells)))
    L += ["\\bottomrule", "\\end{tabular}", "\\end{table*}"]
    write("fd4", L)
# ---------- 10. chi_GaP dependence of r* (round-4 review) ----------
rc = rows(["results", "paper_design_map_chi", "rstar_chi.csv"])
if rc:
    RS = {(float(r["chiGaP"]), r["Condition"], float(r["tauP_ns"])): float(r["r_star"]) for r in rc}
    chis = sorted({k[0] for k in RS}); taus = sorted({k[2] for k in RS})
    fmt = lambda v: "--" if math.isinf(v) else ("%.1f" % v if v >= 10 else "%.2f" % v)
    L = [r"\begin{table*}[htbp]", r"\centering",
         r"\caption{Crossover ratio $r^*(\tau_p)$ of RHJ versus RBSF--GaP as a function of the GaP electron affinity $\chi_{\rm GaP}$ ($\chi_{\rm abs}=4.0$~eV, so $\Delta E_c=4.0~{\rm eV}-\chi_{\rm GaP}$; nominal 3.8~eV), $\mu_n=500$~cm$^2$\,V$^{-1}$s$^{-1}$, $S_{\rm back}=10^5$~cm\,s$^{-1}$, AM1.5G / 6500~K, 1000~lx. ``--'': the RHJ does not reach RBSF--GaP within $r\le60$.}",
         r"\label{tab:supp-rstar-chi}", r"\scriptsize", r"\begin{tabular}{l" + "c" * len(chis) + "}", r"\toprule",
         r"$\chi_{\rm GaP}$ (eV) & " + " & ".join("%.2f" % c if round(c, 1) != c else "%.1f" % c for c in chis) + r" \\",
         r"$\Delta E_c$ (eV) & " + " & ".join(("$%+.2f$" % (4.0 - c)) for c in chis) + r" \\", r"\midrule",
         r"$\tau_p$ (ns) & \multicolumn{%d}{c}{$r^*$, AM1.5G / 1000~lx} \\" % len(chis), r"\midrule"]
    for t in taus:
        L.append("%g & %s \\\\" % (t, " & ".join("%s / %s" % (fmt(RS[(c, "AM1.5G", t)]), fmt(RS[(c, "6500K_1000lx", t)])) for c in chis)))
    L += [r"\bottomrule", r"\end{tabular}", r"\end{table*}"]
    write("chi", L)
print("tables written to", OUT)
