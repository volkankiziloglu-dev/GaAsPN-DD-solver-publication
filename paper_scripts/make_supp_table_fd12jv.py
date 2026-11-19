"""Supplementary table of the F_1/2 (consistent degenerate SG) versus Boltzmann J-V comparison (run_paper_fd12_jv.m).
Input: results/paper_fd12_jv/*_summary.csv. Output: latex/tables/supp_fd12_jv.tex and results/paper_fd12_jv/fd12_jv_summary.csv."""
import glob, os
import pandas as pd

root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = os.path.join(root, "results", "paper_fd12_jv")
files = [f for f in sorted(glob.glob(os.path.join(R, "*_summary.csv"))) if os.path.basename(f) != "fd12_jv_summary.csv"]
d = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
order_c = {"AM1.5G": 0, "6500K_1000lx": 1, "6500K_100lx": 2}
order_a = {"HJ-0": 0, "F-HJ": 1, "RBSF-GaP": 2, "RHJ": 3}
d = d[d.Condition.isin(order_c)].copy()
d["oc"] = d.Condition.map(order_c); d["oa"] = d.Architecture.map(order_a)
d = d.sort_values(["oc", "oa"])
d["dPCE"] = d.PCEFD - d.PCEB
d["dVoc_mV"] = (d.VocFD - d.VocB) * 1e3
d["dJsc_rel"] = d.JscFD / d.JscB - 1
d["dFF_pts"] = (d.FFFD - d.FFB) * 100
d.drop(columns=["oc", "oa"]).to_csv(os.path.join(R, "fd12_jv_summary.csv"), index=False, lineterminator=chr(10))
cname = {"AM1.5G": "AM1.5G", "6500K_1000lx": r"1000~lx", "6500K_100lx": r"100~lx"}
aname = {"HJ-0": "HJ-0", "F-HJ": "F-HJ", "RBSF-GaP": r"RBSF--GaP", "RHJ": "RHJ"}
L = [r"\begin{table}[htbp]", r"\centering",
     r"\caption{Fermi--Dirac ($F_{1/2}$, consistent degenerate Scharfetter--Gummel flux) versus Boltzmann statistics, nominal devices. At every voltage the $F_{1/2}$ solution is started from the Boltzmann solution at the same voltage; both curves are evaluated with the unified extractor on all finite points. Conv.: points meeting the outer-iteration criteria of the respective solver (\texttt{cases/run\_paper\_fd12\_jv.m}, \texttt{results/paper\_fd12\_jv/}).}",
     r"\label{tab:supp-fd12jv}", r"\footnotesize", r"\setlength{\tabcolsep}{4pt}",
     r"\begin{tabular}{llrrrrrr}", r"\hline",
     r"Condition & Architecture & PCE$_{\rm B}$ (\%) & PCE$_{F_{1/2}}$ (\%) & $\Delta$PCE (pts) & $\Delta V_{oc}$ (mV) & $\Delta J_{sc}/J_{sc}$ & Conv.\ B/$F_{1/2}$ \\",
     r"\hline"]
prev = None
for _, r in d.iterrows():
    c = cname[r.Condition] if r.Condition != prev else ""
    if prev is not None and r.Condition != prev:
        L.append(r"\hline")
    prev = r.Condition
    L.append("%s & %s & %.2f & %.2f & %+.3f & %+.1f & %.1e & %d/%d of %d \\\\" % (
        c, aname[r.Architecture], r.PCEB, r.PCEFD, r.dPCE, r.dVoc_mV, r.dJsc_rel, r.nConvB, r.nConvFD, r.nV))
L += [r"\hline", r"\end{tabular}", r"\end{table}"]
out = os.path.join(root, "latex", "tables", "supp_fd12_jv.tex")
open(out, "w", encoding="utf-8", newline=chr(10)).write(chr(10).join(L) + chr(10))
print(d[["Condition", "Architecture", "PCEB", "PCEFD", "dPCE", "dVoc_mV", "nConvB", "nConvFD"]].to_string(index=False))
