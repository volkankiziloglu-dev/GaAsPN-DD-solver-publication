"""[scaps_klitzke_all.py, 2026-11-13: copy of scaps_prod_all.py for the Klitzke v112 devices (sesame_input/*_v112_*).]
SCAPS-1D 3.3.12 cross-check of the PRODUCTION configuration.
Builds (scaps_build_prod.py), runs SCAPS and compares with the MATLAB J-V (all finite points, cases/export_for_scaps_prod.m) on one voltage
grid with one extractor (J_sc at 0 V, linear V_oc, dense-grid P_max), as scaps_all.py does for the SRH-only check.
Usage: python scaps_prod_all.py [tag ...]   Output: scaps_output/scaps_prod_crosscheck.csv
"""
import os, re, sys, json, subprocess
import numpy as np
from scaps_run import run, WORK

here = os.path.dirname(os.path.abspath(__file__))
TAGS = ["KHJ_v112_100lx", "KHJ_v112_1000lx", "KHJrad_v112_100lx", "KHJrad_v112_1000lx", "KRHJlo_v112_100lx", "KRHJlo_v112_1000lx"]
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
        if m: rows.append((float(m.group(1)), float(m.group(2)) * 1e-3))
    return np.array(rows)


out = []
for tag in tags:
    subprocess.run([sys.executable, "-I", os.path.join(here, "scaps_build_klitzke.py"), tag], check=True, stdout=subprocess.DEVNULL,
                   env=dict(os.environ))
    if True:
        p = os.path.join(WORK, tag + ".def"); t = open(p, encoding="latin-1", newline="").read()
        t = re.sub(r"^generation limit :.*", "generation limit :  1.00e-06 microA/cm2", t, flags=re.M)
        t = re.sub(r"^recombination limit :.*", "recombination limit :  1.00e-09 microA/cm2", t, flags=re.M)
        open(p, "w", encoding="latin-1", newline="").write(t)
    meta = json.load(open(os.path.join(here, "sesame_input", tag + "_meta.json")))
    gen = os.path.join(WORK, tag + ".gen")
    rc, log, ivp, dt = run(tag + "_scaps", ["load generationfile " + gen, "action generationfromfile", "action light"],
                           defname=tag, vstart=0.0, vstop=1.60, dv=0.02, stopvoc=0, timeout=900)
    if not ivp:
        print(tag, "NO IV", rc, log[:300]); continue
    sc = read_iv(ivp)
    if len(sc) < 5:
        print(tag, "too few points", len(sc), log[:300]); continue
    ml = np.genfromtxt(os.path.join(here, "sesame_input", tag + "_mlJV.csv"), delimiter=",", names=True)
    ok = (ml["Accepted"] == 1) & np.isfinite(ml["J_Acm2"])
    V = np.arange(0, 1.5601, 0.02)   # the 0.02 V grid of both codes (MPP lies below 1.3 V at 100-200 lx)
    V = V[(V <= sc[-1, 0] + 1e-9) & (V <= ml["V_V"][ok].max() + 1e-9)]
    Js = np.interp(V, sc[:, 0], sc[:, 1]); Jm = np.interp(V, ml["V_V"][ok], ml["J_Acm2"][ok])
    Pin = meta["Pin_Wm2"]
    try:
        mS = metrics(V, Js, Pin); mM = metrics(V, Jm, Pin)
    except Exception as e:
        print(tag, "metrics failed", e, "last V", sc[-1, 0]); continue
    out.append((tag, mS, mM))
    print("%-22s SCAPS Jsc=%.5g Voc=%.2f FF=%.3f PCE=%.4f | MATLAB Jsc=%.5g Voc=%.2f FF=%.3f PCE=%.4f  (%.0fs)" % (
        (tag,) + (mS[0] * 1e3, mS[1] * 1e3, mS[2] * 100, mS[3]) + (mM[0] * 1e3, mM[1] * 1e3, mM[2] * 100, mM[3]) + (dt,)))

with open(os.path.join(here, "scaps_output", "scaps_klitzke_v112_crosscheck.csv"), "w") as f:
    f.write("tag,code,Jsc_mAcm2,Voc_mV,FF_pct,PCE_pct\n")
    for tag, mS, mM in out:
        for nm, m in (("SCAPS", mS), ("MATLAB", mM)):
            f.write("%s,%s,%.5f,%.2f,%.3f,%.4f\n" % (tag, nm, m[0] * 1e3, m[1] * 1e3, m[2] * 100, m[3]))
