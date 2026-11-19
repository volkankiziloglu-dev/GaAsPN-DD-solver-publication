"""Summarise the Sesame/MATLAB Klitzke-type cross-check (run after sesame_klitzke_crosscheck.py): writes sesame_output/klitzke_crosscheck_summary.csv and prints LaTeX rows."""
import numpy as np, json, os
here = os.path.dirname(os.path.abspath(__file__))
def metrics(Vv, Jv, Pin):
    Jsc = np.interp(0, Vv, Jv); idx = np.where((Jv[:-1] > 0) & (Jv[1:] <= 0))[0]; i = idx[0]
    Voc = Vv[i] + (0 - Jv[i]) * (Vv[i+1] - Vv[i]) / (Jv[i+1] - Jv[i])
    Vf = np.linspace(0, Voc, 4000); P = (Vf * np.interp(Vf, Vv, Jv)).max()
    return Jsc, Voc, P / (Jsc * Voc), 100 * P * 1e4 / Pin
names = {"KHJ": "HJ, SRH", "KHJrad": "HJ, SRH+rad+Auger", "KRHJhi": r"RHJ, $N_D=5\times10^{17}$", "KRHJlo": r"RHJ, $N_D=2\times10^{16}$"}
rows = []; out = ["tag,dJsc_pct,dVoc_mV,dFF_pp,dPCE_pp,maxrelJ_below_Voc"]
for dv in names:
    for lx in (100, 1000):
        tag = f"{dv}_{lx}lx"
        meta = json.load(open(f"sesame_input/{tag}_meta.json"))
        ml = np.genfromtxt(f"sesame_input/{tag}_mlJV.csv", delimiter=",", names=True); ok = ml["Accepted"] == 1
        se = np.genfromtxt(f"sesame_output/{tag}_sesame_JV.csv", delimiter=",", names=True)
        V = se["V"]; Js = se["J_Acm2"]; Jm = ml["J_Acm2"][ok]
        s = metrics(V, Js, meta["Pin_Wm2"]); m = metrics(V, Jm, meta["Pin_Wm2"])
        below = V <= m[1] + 1e-9
        rel = np.max(np.abs(Js[below] - Jm[below]) / m[0])
        d = (100*(s[0]/m[0]-1), (s[1]-m[1])*1e3, (s[2]-m[2])*100, s[3]-m[3])
        out.append("%s,%.2e,%.4f,%.4f,%.4f,%.2e" % ((tag,) + d + (rel,)))
        rows.append((names[dv], lx, s, m, d, rel))
        print("%-24s %4d lx  Jsc %.3f/%.3f uA  Voc %.1f/%.1f  FF %.2f/%.2f  PCE %.3f/%.3f   dJsc=%.1e%% dVoc=%.3f mV dFF=%.4f dPCE=%.4f  max|dJ|/Jsc(V<=Voc)=%.1e" % (
            names[dv], lx, s[0]*1e6, m[0]*1e6, s[1]*1e3, m[1]*1e3, s[2]*100, m[2]*100, s[3], m[3], d[0], d[1], d[2], d[3], rel))
open("sesame_output/klitzke_crosscheck_summary.csv", "w").write("\n".join(out) + "\n")
