"""Summary of the Sesame/MATLAB and SCAPS/MATLAB cross-check of the Klitzke v112 devices (cases/export_klitzke_for_sesame_v112.m;
sesame_klitzke_crosscheck.py; scaps_klitzke_all.py). Writes sesame_output/klitzke_v112_crosscheck_summary.csv and prints LaTeX rows."""
import os, json
import numpy as np
here = os.path.dirname(os.path.abspath(__file__))


def met(V, J, Pin):
    Jsc = np.interp(0, V, J); i = np.where((J[:-1] > 0) & (J[1:] <= 0))[0][0]
    Voc = V[i] - J[i] * (V[i + 1] - V[i]) / (J[i + 1] - J[i])
    Vf = np.linspace(0, Voc, 4000); P = (Vf * np.interp(Vf, V, J)).max()
    return Jsc, Voc, P / (Jsc * Voc), 100 * P * 1e4 / Pin


nm = {"KHJ": "HJ, SRH only", "KHJrad": "HJ, SRH + radiative + Auger", "KRHJlo": r"RHJ, $N_D=2\times10^{16}$~cm$^{-3}$"}
out = ["tag,dJsc_pct,dVoc_mV,dFF_pp,dPCE_pp,maxrelJ"]
for d in nm:
    for lx in (100, 1000):
        t = "%s_v112_%dlx" % (d, lx); me = json.load(open(os.path.join(here, "sesame_input", t + "_meta.json")))
        ml = np.genfromtxt(os.path.join(here, "sesame_input", t + "_mlJV.csv"), delimiter=",", names=True); ok = ml["Accepted"] == 1
        se = np.genfromtxt(os.path.join(here, "sesame_output", t + "_sesame_JV.csv"), delimiter=",", names=True)
        V = se["V"]; Js = se["J_Acm2"]; Jm = np.interp(V, ml["V_V"][ok], ml["J_Acm2"][ok])
        s = met(V, Js, me["Pin_Wm2"]); m = met(V, Jm, me["Pin_Wm2"]); b = V <= m[1] + 1e-9
        r = np.max(np.abs(Js[b] - Jm[b])) / m[0]
        out.append("%s,%.3e,%.4f,%.4f,%.4f,%.2e" % (t, 100 * (s[0] / m[0] - 1), (s[1] - m[1]) * 1e3, (s[2] - m[2]) * 100, s[3] - m[3], r))
        e = int(np.floor(np.log10(r))); mant = int(round(r / 10 ** e))
        print("%s & %d & %.3f/%.3f & %.1f/%.1f & %.2f/%.2f & %.3f/%.3f & $%d\\times10^{%d}$ \\\\" % (
            nm[d], lx, s[0] * 1e6, m[0] * 1e6, s[1] * 1e3, m[1] * 1e3, s[2] * 100, m[2] * 100, s[3], m[3], mant, e))
open(os.path.join(here, "sesame_output", "klitzke_v112_crosscheck_summary.csv"), "w").write("\n".join(out) + "\n")
