"""Independent Sesame (NIST) cross-check of the paper architectures RBSF-GaP and RHJ under AM1.5G.

Inputs come from cases/export_for_sesame.m: per-node material parameters, doping and the MATLAB generation profile (identical optical input),
SRH only, ideal heterointerfaces (frontIRV = rearIRV = 0), rear minority-carrier velocity Sback. The same x nodes are used as the MATLAB mesh.
Run:  python sesame_hetero_crosscheck.py RBSF-GaP   (or RHJ)
Output: external_validation/sesame_output/<code>_sesame_JV.csv and a printed comparison with the MATLAB J-V.
"""
import sys, os, json
import numpy as np
import sesame

code = sys.argv[1] if len(sys.argv) > 1 else "RBSF-GaP"
here = os.path.dirname(os.path.abspath(__file__))
inp = os.path.join(here, "sesame_input")
outd = os.path.join(here, "sesame_output"); os.makedirs(outd, exist_ok=True)

nodes = np.genfromtxt(os.path.join(inp, f"{code}_nodes.csv"), delimiter=",", names=True)
meta = json.load(open(os.path.join(inp, f"{code}_meta.json")))
ml = np.genfromtxt(os.path.join(inp, f"{code}_mlJV.csv"), delimiter=",", names=True)

x = nodes["x_cm"]
layers = np.unique(nodes["layer"]).astype(int)
# boundaries between layers: midpoint between last node of a layer and first node of the next
bounds = {}
for k, kn in zip(layers[:-1], layers[1:]):
    xl = x[nodes["layer"] == k].max(); xf = x[nodes["layer"] == kn].min()
    bounds[(k, kn)] = 0.5 * (xl + xf)

def loc(k):
    lo = -np.inf if k == layers[0] else bounds[(layers[list(layers).index(k) - 1], k)]
    hi = np.inf if k == layers[-1] else bounds[(k, layers[list(layers).index(k) + 1])]
    return lambda pos: (pos > lo) & (pos <= hi) if lo > -np.inf else (pos <= hi)

system = sesame.Builder(x)
for k in layers:
    r = nodes[nodes["layer"] == k][0]
    mat = {"Nc": r["Nc"], "Nv": r["Nv"], "Eg": r["Eg"], "affinity": r["chi"], "epsilon": r["eps"],
           "mu_e": r["muN"], "mu_h": r["muP"], "tau_e": r["tau"], "tau_h": r["tau"], "Et": 0, "B": 0}
    system.add_material(mat, loc(k))
    if r["Nd"] > 0:
        system.add_donor(r["Nd"], loc(k))
    if r["Na"] > 0:
        system.add_acceptor(r["Na"], loc(k))

system.contact_type("Ohmic", "Ohmic")
# front: ideal for both carriers (n+ GaP window); rear (p+ GaP): majority ideal, minority (electron) velocity Sback
system.contact_S(1e7, 1e7, meta["Sback"], 1e7)
G = nodes["G"]
system.generation(lambda xx: np.interp(xx, x, G))

Nn = x.size
poisson_eq = sesame.solve(system, compute="Poisson", guess=None, htp=5, verbose=False)
guess0 = {"efn": np.zeros(Nn), "efp": np.zeros(Nn), "v": poisson_eq["v"]}
eq = sesame.solve(system, compute="all", guess=guess0, htp=5, verbose=False)

V = np.concatenate([np.arange(0, 1.30, 0.1), np.arange(1.30, 1.56, 0.02)])
J = sesame.IVcurve(system, V, os.path.join(outd, f"{code}_sesame"), guess=eq, fmt="npz", verbose=False, htp=3)
J = np.asarray(J) * system.scaling.current          # A/cm2
sgn = 1.0 if J[0] > 0 else -1.0
J = sgn * J
np.savetxt(os.path.join(outd, f"{code}_sesame_JV.csv"), np.c_[V, J], delimiter=",", header="V,J_Acm2", comments="")

mlV = ml["V_V"][ml["Accepted"] == 1]; mlJ = ml["J_Acm2"][ml["Accepted"] == 1]
Jml = np.interp(V, mlV, mlJ)
print(f"== {code}: Sesame vs MATLAB (AM1.5G, SRH only, ideal interfaces)")
print(" V(V)   J_Sesame(mA/cm2)  J_MATLAB(mA/cm2)  diff(mA/cm2)")
for v, a, b in zip(V, J, Jml):
    print(f"{v:5.2f}  {a*1e3:14.4f}  {b*1e3:14.4f}  {(a-b)*1e3:+10.4f}")

def metrics(Vv, Jv, Pin):
    Jsc = np.interp(0, Vv, Jv)
    idx = np.where((Jv[:-1] > 0) & (Jv[1:] <= 0))[0]
    if len(idx) == 0:
        return Jsc, np.nan, np.nan, np.nan
    i = idx[0]; Voc = Vv[i] + (0 - Jv[i]) * (Vv[i+1] - Vv[i]) / (Jv[i+1] - Jv[i])
    Vf = np.linspace(0, Voc, 4000); Jf = np.interp(Vf, Vv, Jv); P = (Vf * Jf).max()
    return Jsc, Voc, P / (Jsc * Voc), 100 * P * 1e4 / Pin
Pin = meta["Pin_Wm2"]
s = metrics(V, J, Pin); m = metrics(V, Jml, Pin)
print(f"Sesame: Jsc={s[0]*1e3:.3f} mA/cm2 Voc={s[1]*1e3:.1f} mV FF={s[2]*100:.2f}% PCE={s[3]:.3f}%  (grid-limited near MPP)")
print(f"MATLAB: Jsc={m[0]*1e3:.3f} mA/cm2 Voc={m[1]*1e3:.1f} mV FF={m[2]*100:.2f}% PCE={m[3]:.3f}%  (same voltage grid)")
print(f"MATLAB full-sweep metrics: Jsc={meta['Jsc_mAcm2']:.3f} Voc={meta['Voc_mV']:.1f} FF={meta['FF_pct']:.2f} PCE={meta['PCE_pct']:.3f}")
