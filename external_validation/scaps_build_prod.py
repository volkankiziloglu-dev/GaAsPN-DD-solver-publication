"""SCAPS-1D 3.3.12 input for the PRODUCTION configuration of a paper architecture (simulated referee R1-1, docs/HAKEM_SIMULASYONU_2026-11-06.md).

Same construction as scaps_build.py (layers, flat-band contacts, midgap SRH, MATLAB generation file), plus the production recombination:
  * radiative K_rad = B_rad and Auger c_n = c_p (per layer, from the exported nodes),
  * an interface defect at every heterointerface with the MATLAB interface recombination velocity S_int: neutral, narrow uniform band
    (width 0.01 eV) at the middle of the interface gap, sigma_n = sigma_p = 1e-19 m^2 on both sides, v_th = 1e5 m/s, so that
    sigma * v_th * N_s = S_int (N_s = 8.91e15 m^-2 for 8910 cm/s).
SCAPS evaluates the interface recombination with its own interface densities; MATLAB uses the geometric mean of the densities of the
two nodes adjacent to the interface, so the comparison tests that formulation.
Input : scaps_prod_input/<tag>_nodes.csv, <tag>_meta.json (cases/export_for_scaps_prod.m)
Output: <SCAPS_WORK>/<tag>.def, <tag>.gen, <tag>.script
Environment: SCAPS_WORK (default ~/scaps_xcheck), SCAPS_TEMPLATE (definition file used as the text template).
"""
import os, re, sys, json
import numpy as np

tag = sys.argv[1] if len(sys.argv) > 1 else "P_RBSF-GaP"
here = os.path.dirname(os.path.abspath(__file__))
work = os.environ.get("SCAPS_WORK", os.path.join(os.path.expanduser("~"), "scaps_xcheck")); os.makedirs(work, exist_ok=True)
tmpl_path = os.environ["SCAPS_TEMPLATE"]   # SCAPS-1D definition file used as the text template (not distributed; see docs/data_sources.md)
T = open(tmpl_path, encoding="latin-1").read().split("\n")
T = [l.rstrip("\r") for l in T]

nodes = np.genfromtxt(os.path.join(here, "scaps_prod_input", f"{tag}_nodes.csv"), delimiter=",", names=True)
meta = json.load(open(os.path.join(here, "scaps_prod_input", f"{tag}_meta.json")))
x = nodes["x_cm"] * 1e-2
layers = list(dict.fromkeys(nodes["layer"].astype(int)))
layer_info = []
for k in layers:
    sel = nodes["layer"] == k
    xs = x[sel]
    layer_info.append((k, int(sel.sum()), float(np.median(np.diff(xs))), nodes[sel][0]))
if len(layers) > 6:
    raise SystemExit("this builder handles up to 6 layers")
Sint = meta["Sint_cmps"]; Sint = [Sint] if np.isscalar(Sint) else list(Sint)


def find_block(lines, start_pat):
    for i, l in enumerate(lines):
        if l.startswith(start_pat):
            return i
    raise KeyError(start_pat)


i_back = find_block(T, "back contact"); i_layer1 = find_block(T, "layer")
i_front = find_block(T, "front contact"); i_work = find_block(T, "working point")
header = T[:i_back]
idx_layers = [i for i, l in enumerate(T) if l == "layer"]
layer_tmpl = T[idx_layers[-1]:i_front]
back_tmpl = T[i_back:i_layer1]; front_tmpl = T[i_front:i_work]; work_tmpl = T[i_work:]
num_re = re.compile(r"[-+]?\d+\.?\d*(?:[eE][-+]?\d+)?")


def set_vals(lines, key, value, fmt="%.6e"):
    for i, l in enumerate(lines):
        if l.startswith(key + " :") or l.startswith(key + ":"):
            head, tail = l.split(":", 1)
            toks = list(num_re.finditer(tail)); new = tail
            for j in sorted([0, 5, 6], reverse=True):
                if j < len(toks):
                    m = toks[j]; new = new[:m.start()] + (fmt % value) + new[m.end():]
            lines[i] = head + ":" + new
            return
    raise KeyError(key)


def set_scalar(lines, key, value, fmt="%.6e"):
    for i, l in enumerate(lines):
        if l.startswith(key):
            head, tail = l.split(":", 1)
            m = num_re.search(tail)
            lines[i] = head + ":" + tail[:m.start()] + (fmt % value) + tail[m.end():]
            return
    raise KeyError(key)


def make_layer(name, d_m, r):
    L = list(layer_tmpl)
    L[2] = "name : " + name
    set_scalar(L, "d :", d_m)
    set_vals(L, "eps", r["eps"], "%.6f"); set_vals(L, "chi", r["chi"], "%.6f"); set_vals(L, "Eg", r["Eg"], "%.6f")
    set_vals(L, "Nc", r["Nc"] * 1e6); set_vals(L, "Nv", r["Nv"] * 1e6)
    set_vals(L, "mu_n", r["muN"] * 1e-4); set_vals(L, "mu_p", r["muP"] * 1e-4)
    set_vals(L, "K_rad", r["Brad"] * 1e-6)                 # cm3/s -> m3/s
    set_vals(L, "c_n_auger", r["Cn"] * 1e-12); set_vals(L, "c_p_auger", r["Cp"] * 1e-12)   # cm6/s -> m6/s
    set_vals(L, "Na(uniform)", r["Na"] * 1e6); set_vals(L, "Nd(uniform)", r["Nd"] * 1e6)
    sigma = 1e-19; vth = 1e5; Nt = 1.0 / (r["tau"] * sigma * vth)
    set_vals(L, "Nt(uniform)", Nt); set_scalar(L, "Nleft", Nt); set_scalar(L, "Nright", Nt)
    for i, l in enumerate(L):
        if l.startswith("sigma_n"): L[i] = "sigma_n : 1.000e-19\t[m^2]"
        if l.startswith("sigma_p"): L[i] = "sigma_p : 1.000e-19\t[m^2]"
        if l.startswith("Reference for defect energy"):
            L[i] = "Reference for defect energy :  1 ; 0 = above EV (SCAPS<2.7); 1 = above Ei; 2 = below EC"
        if l.startswith("Et :"): L[i] = "Et :   0.000\t[eV]"
        if l.startswith("Ekar"): L[i] = "Ekar :  0.0000 [eV]"
    return L


def interface_block(name, S_cmps):
    out = ["interface properties", "", "interfacename : " + name, "Intraband tunneling : 0 (0: no tunneling, 1: with tunneling)",
           "Relative electron mass :  1.00e+00", "Relative hole mass :  1.00e+00", ""]
    if S_cmps > 0:
        width = 0.01; Ns = (S_cmps * 1e-2) / (1e-19 * 1e5)          # m^-2
        out += ["interface recombination", "", "type : neutral",
                "sigma_nleft : 1.000e-19 [m^2]", "sigma_nright : 1.000e-19 [m^2]",
                "sigma_pleft : 1.000e-19 [m^2]", "sigma_pright : 1.000e-19 [m^2]",
                "Tunneling to trap: 0 (0: no tunneling, 1: with tunneling)",
                "Relative electron mass :  1.00e+00", "Relative hole mass :  1.00e+00",
                "energy distribution : uniform",
                "Reference for defect energy :  9 ; 3 = above highest EV; 4 = below lowest EC; 5 = above Ei left; 6 = above Ei right; 7 = above EV left (= SCAPS%<2.7 reference); 8 = above EV right; 9 = above middle of interface gap",
                "Et :  0.00 [eV]", "Ekar :  %.2f [eV]" % width, "N : %.3e [/m^2.eV]" % (Ns / width), ""]
    out += [""]
    return out


kT_eV = 0.025852
def flat_phi(r):
    Nd_, Na_ = r["Nd"], r["Na"]
    if Nd_ >= Na_:
        return r["chi"] + kT_eV * np.log(r["Nc"] / max(Nd_ - Na_, 1.0))
    return r["chi"] + r["Eg"] - kT_eV * np.log(r["Nv"] / (Na_ - Nd_))


out = ["> SCAPS 3.3.12 production-configuration cross-check of %s (exported from the MATLAB model)." % tag,
       "> Midgap SRH, radiative, Auger and interface defects with S_int; generation from file (MATLAB G(x))."]
out += [l for l in header if not l.startswith(">")]
bl = list(back_tmpl)
set_scalar(bl, "Sn", 1e5); set_scalar(bl, "Sp", 1e5); set_scalar(bl, "Fi_m", flat_phi(layer_info[0][3]), "%.6f")
for i, l in enumerate(bl):
    if l.startswith("flatband"): bl[i] = "flatband : 1 (1: flatband contact imposed; 0: fixed Phi_m value)"
out += bl
for n, (k, cnt, dx, r) in enumerate(layer_info):
    out += make_layer("L%d" % k, cnt * dx, r)
    if n < len(layer_info) - 1:
        S = Sint[n] if n < len(Sint) else 0.0
        same = abs(r["chi"] - layer_info[n + 1][3]["chi"]) < 1e-9 and abs(r["Eg"] - layer_info[n + 1][3]["Eg"]) < 1e-9
        out += interface_block("L%d / L%d" % (k, layer_info[n + 1][0]), 0.0 if same else S)
fl = list(front_tmpl)
set_scalar(fl, "Sn", meta["Sback"] * 1e-2); set_scalar(fl, "Sp", 1e5); set_scalar(fl, "Fi_m", flat_phi(layer_info[-1][3]), "%.6f")
for i, l in enumerate(fl):
    if l.startswith("flatband"): fl[i] = "flatband : 1 (1: flatband contact imposed; 0: fixed Phi_m value)"
out += fl + work_tmpl
open(os.path.join(work, f"{tag}.def"), "w", encoding="latin-1").write("\n".join(out))

dxn = np.empty_like(x)
for k, cnt, dx, r in layer_info:
    dxn[nodes["layer"] == k] = dx
xcen = x + 0.5 * dxn
g = nodes["G"] * 1e6
Ltot = sum(cnt * dx for k, cnt, dx, r in layer_info)
EPS = 1e-3
pts = [(0.0, g[0])]
pts_layers = [np.where(nodes["layer"] == k)[0] for k in layers]
for n_, idx in enumerate(pts_layers):
    pts += [(xcen[i] * 1e6, g[i]) for i in idx]
    if n_ < len(pts_layers) - 1:
        bnd = x[pts_layers[n_ + 1][0]] * 1e6
        pts += [(bnd - EPS, g[idx[-1]]), (bnd + EPS, g[pts_layers[n_ + 1][0]])]
pts.append((Ltot * 1e6, g[-1])); pts.sort(key=lambda p: p[0])
with open(os.path.join(work, f"{tag}.gen"), "w", newline="") as f:
    f.write("Generation profile exported from MATLAB\r\n\r\nx (um)\tGeh (#/m3.s)\r\n")
    for a, b in pts:
        f.write("%.8f\t%.5E\r\n" % (a, b))
print("built", tag, "layers:", [(int(k), cnt, dx) for k, cnt, dx, r in layer_info], "Sint:", Sint)
