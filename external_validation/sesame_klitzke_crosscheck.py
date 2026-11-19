"""Independent Sesame (NIST) cross-check of Klitzke-type GaInP homojunction / rear-heterojunction devices under the CW-LED spectrum.

Inputs from cases/export_klitzke_for_sesame.m: per-node parameters incl. separate electron/hole SRH lifetimes (one centre at E_i), radiative B and Auger Cn, Cp,
and the MATLAB generation profile. Ideal heterointerfaces, ideal front contact, minority velocity Sback at the rear.
Run:  python sesame_klitzke_crosscheck.py KHJ_1000lx [more tags]   -> sesame_output/<tag>_sesame_JV.csv and a metrics line per tag.
"""
import sys, os, json
import numpy as np
import sesame

here = os.path.dirname(os.path.abspath(__file__))
inp = os.path.join(here, "sesame_input"); outd = os.path.join(here, "sesame_output"); os.makedirs(outd, exist_ok=True)


def metrics(Vv, Jv, Pin):
    Jsc = np.interp(0, Vv, Jv)
    idx = np.where((Jv[:-1] > 0) & (Jv[1:] <= 0))[0]
    if len(idx) == 0:
        return Jsc, np.nan, np.nan, np.nan
    i = idx[0]; Voc = Vv[i] + (0 - Jv[i]) * (Vv[i + 1] - Vv[i]) / (Jv[i + 1] - Jv[i])
    Vf = np.linspace(0, Voc, 4000); Jf = np.interp(Vf, Vv, Jv); P = (Vf * Jf).max()
    return Jsc, Voc, P / (Jsc * Voc), 100 * P * 1e4 / Pin


def run(tag):
    nodes = np.genfromtxt(os.path.join(inp, f"{tag}_nodes.csv"), delimiter=",", names=True)
    meta = json.load(open(os.path.join(inp, f"{tag}_meta.json")))
    ml = np.genfromtxt(os.path.join(inp, f"{tag}_mlJV.csv"), delimiter=",", names=True)
    x = nodes["x_cm"]; layers = np.unique(nodes["layer"]).astype(int)
    bounds = {}
    for k, kn in zip(layers[:-1], layers[1:]):
        bounds[(k, kn)] = 0.5 * (x[nodes["layer"] == k].max() + x[nodes["layer"] == kn].min())

    def loc(k):
        i = list(layers).index(k)
        lo = -np.inf if i == 0 else bounds[(layers[i - 1], k)]
        hi = np.inf if i == len(layers) - 1 else bounds[(k, layers[i + 1])]
        return lambda pos: (pos > lo) & (pos <= hi) if lo > -np.inf else (pos <= hi)

    system = sesame.Builder(x)
    for k in layers:
        r = nodes[nodes["layer"] == k][0]
        mat = {"Nc": r["Nc"], "Nv": r["Nv"], "Eg": r["Eg"], "affinity": r["chi"], "epsilon": r["eps"],
               "mu_e": r["muN"], "mu_h": r["muP"], "tau_e": r["tauN"], "tau_h": r["tauP"], "Et": 0,
               "B": r["Brad"], "Cn": r["Cn"], "Cp": r["Cp"]}
        system.add_material(mat, loc(k))
        if r["Nd"] > 0: system.add_donor(r["Nd"], loc(k))
        if r["Na"] > 0: system.add_acceptor(r["Na"], loc(k))
    system.contact_type("Ohmic", "Ohmic")
    system.contact_S(1e7, 1e7, meta["Sback"], 1e7)
    G = nodes["G"]
    system.generation(lambda xx: np.interp(xx, x, G))
    Nn = x.size
    poisson_eq = sesame.solve(system, compute="Poisson", guess=None, htp=5, verbose=False)
    guess0 = {"efn": np.zeros(Nn), "efp": np.zeros(Nn), "v": poisson_eq["v"]}
    eq = sesame.solve(system, compute="all", guess=guess0, htp=5, verbose=False)
    V = ml["V_V"][ml["Accepted"] == 1]
    J = sesame.IVcurve(system, V, os.path.join(outd, f"{tag}_sesame"), guess=eq, fmt="npz", verbose=False, htp=3)
    J = np.asarray(J) * system.scaling.current
    J = (1.0 if J[0] > 0 else -1.0) * J
    np.savetxt(os.path.join(outd, f"{tag}_sesame_JV.csv"), np.c_[V, J], delimiter=",", header="V,J_Acm2", comments="")
    Jml = ml["J_Acm2"][ml["Accepted"] == 1]
    Pin = meta["Pin_Wm2"]
    s = metrics(V, J, Pin); m = metrics(V, Jml, Pin)
    print(f"{tag}: Sesame Jsc={s[0]*1e6:.4f} uA/cm2 Voc={s[1]*1e3:.1f} mV FF={s[2]*100:.2f} PCE={s[3]:.3f} | "
          f"MATLAB Jsc={m[0]*1e6:.4f} Voc={m[1]*1e3:.1f} FF={m[2]*100:.2f} PCE={m[3]:.3f}", flush=True)


for t in sys.argv[1:]:
    try:
        run(t)
    except Exception as e:
        print(t, "FAILED", repr(e)[:200], flush=True)
