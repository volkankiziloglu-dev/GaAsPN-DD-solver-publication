"""SCAPS-1D 3.3.12 cross-check of the paper architectures: build, run, compare with MATLAB and Sesame on the same voltage grid.
Usage: python scaps_all.py [tag ...]   (default: all tags below). Output: scaps_output/scaps_crosscheck.csv
"""
import os, re, sys, json, subprocess
import numpy as np
from scaps_run import run, WORK

here = os.path.dirname(os.path.abspath(__file__))
TAGS = ["RBSF-GaP", "RHJ", "RBSF-GaP_tau10ns", "RHJ_tau10ns", "RHJ_tau30ns", "RBSF-GaP_indoor", "RHJ_indoor", "F-HJ", "F-HJ_indoor"]
tags = sys.argv[1:] or TAGS
os.makedirs(os.path.join(here, "scaps_output"), exist_ok=True)


def metrics(Vv, Jv, Pin):
    Jsc = np.interp(0, Vv, Jv)
    idx = np.where((Jv[:-1] > 0) & (Jv[1:] <= 0))[0]
    i = idx[0]
    Voc = Vv[i] + (0 - Jv[i]) * (Vv[i + 1] - Vv[i]) / (Jv[i + 1] - Jv[i])
    Vf = np.linspace(0, Voc, 4000); Jf = np.interp(Vf, Vv, Jv); P = (Vf * Jf).max()
    return Jsc, Voc, P / (Jsc * Voc), 100 * P * 1e4 / Pin


def read_iv(path):
    rows = []
    for l in open(path, encoding="latin-1"):
        m = re.match(r"^\s+(-?\d\.\d+)\s+(-?\d\.\d+e[+-]\d+)", l)
        if m: rows.append((float(m.group(1)), float(m.group(2)) * 1e-3))   # J in A/cm2
    return np.array(rows)


out = []
for tag in tags:
    subprocess.run([sys.executable, os.path.join(here, "scaps_build.py"), tag], check=True, stdout=subprocess.DEVNULL)
    indoor = "indoor" in tag
    if indoor:   # tighten SCAPS' absolute current convergence limits for the ~0.1 mA/cm2 indoor currents
        p = os.path.join(WORK, tag + ".def"); t = open(p, encoding="latin-1", newline="").read()
        t = re.sub(r"^generation limit :.*", "generation limit :  1.00e-06 microA/cm2", t, flags=re.M)
        t = re.sub(r"^recombination limit :.*", "recombination limit :  1.00e-09 microA/cm2", t, flags=re.M)
        open(p, "w", encoding="latin-1", newline="").write(t)
    meta = json.load(open(os.path.join(here, "sesame_input", tag + "_meta.json")))
    gen = os.path.join(WORK, tag + ".gen")
    rc, log, ivp, dt = run(tag + "_scaps", ["load generationfile " + gen, "action generationfromfile", "action light"],
                           defname=tag, vstart=0.0, vstop=1.56, dv=0.02, stopvoc=0, timeout=900)
    if not ivp:
        print(tag, "NO IV", log[:200]); continue
    sc = read_iv(ivp)
    if len(sc) < 5:
        print(tag, "too few points", len(sc), log[:200]); continue
    # same voltage grid as the Sesame cross-check
    V = np.concatenate([np.arange(0, 1.30, 0.1), np.arange(1.30, 1.56, 0.02)])
    V = V[V <= sc[-1, 0] + 1e-9]
    Js = np.interp(V, sc[:, 0], sc[:, 1])
    ses = np.genfromtxt(os.path.join(here, "sesame_output", tag + "_sesame_JV.csv"), delimiter=",", names=True)
    ml = np.genfromtxt(os.path.join(here, "sesame_input", tag + "_mlJV.csv"), delimiter=",", names=True)
    ok = ml["Accepted"] == 1
    Jm = np.interp(V, ml["V_V"][ok], ml["J_Acm2"][ok])
    Jse = np.interp(V, ses["V"], ses["J_Acm2"])
    Pin = meta["Pin_Wm2"]
    try:
        mS = metrics(V, Js, Pin)
    except Exception as e:
        print(tag, "SCAPS metrics failed (no Voc reached?)", e, "last V", sc[-1, 0]); continue
    mM = metrics(V, Jm, Pin); mE = metrics(V, Jse, Pin)
    out.append((tag, mS, mE, mM, sc[-1, 0]))
    print("%-18s SCAPS Jsc=%.4g Voc=%.1f FF=%.2f PCE=%.3f | Sesame Jsc=%.4g Voc=%.1f FF=%.2f PCE=%.3f | MATLAB Jsc=%.4g Voc=%.1f FF=%.2f PCE=%.3f  (last V %.2f, %.0fs)" % (
        (tag,) + (mS[0] * 1e3, mS[1] * 1e3, mS[2] * 100, mS[3]) + (mE[0] * 1e3, mE[1] * 1e3, mE[2] * 100, mE[3]) + (mM[0] * 1e3, mM[1] * 1e3, mM[2] * 100, mM[3]) + (sc[-1, 0], dt)))

with open(os.path.join(here, "scaps_output", "scaps_crosscheck.csv"), "w") as f:
    f.write("tag,code,Jsc_mAcm2,Voc_mV,FF_pct,PCE_pct\n")
    for tag, mS, mE, mM, _ in out:
        for nm, m in (("SCAPS", mS), ("Sesame", mE), ("MATLAB", mM)):
            f.write("%s,%s,%.5f,%.2f,%.3f,%.4f\n" % (tag, nm, m[0] * 1e3, m[1] * 1e3, m[2] * 100, m[3]))
