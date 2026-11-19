"""Build SCAPS-1D 3.3.12 input (definition, generation, script) for a paper architecture from the node data exported by cases/export_for_sesame.m.

Same physics as the MATLAB/Sesame cross-check: SRH only (midgap, equal capture lifetimes), ideal heterointerfaces (no interface defects), AM1.5G generation
supplied as the MATLAB generation profile G(x) through a SCAPS generation file, front contact ideal, rear contact with minority-carrier velocity Sback.
Usage: python scaps_build.py <tag>    -> external_validation/scaps_work/<tag>.def, <tag>.gen, <tag>.script
"""
import os, re, sys, json
import numpy as np

tag = sys.argv[1] if len(sys.argv) > 1 else "RBSF-GaP"
here = os.path.dirname(os.path.abspath(__file__))
work = os.environ.get("SCAPS_WORK", os.path.join(os.path.expanduser("~"), "scaps_xcheck")); os.makedirs(work, exist_ok=True)
tmpl_path = os.environ["SCAPS_TEMPLATE"]   # SCAPS-1D definition file used as the text template (not distributed; see docs/data_sources.md)
T = open(tmpl_path, encoding="latin-1").read().split("\n")

nodes = np.genfromtxt(os.path.join(here, "sesame_input", f"{tag}_nodes.csv"), delimiter=",", names=True)
meta = json.load(open(os.path.join(here, "sesame_input", f"{tag}_meta.json")))
x = nodes["x_cm"] * 1e-2          # m
layers = np.unique(nodes["layer"]).astype(int)
dx_all = np.diff(nodes["x_cm"])    # cm
layer_info = []
for k in layers:
    sel = nodes["layer"] == k
    xs = x[sel]; n = sel.sum()
    dx = np.median(np.diff(xs)) if n > 1 else None
    layer_info.append((k, n, dx, nodes[sel][0]))
if len(layers) > 6:
    raise SystemExit("this builder handles up to 6 layers (use the collapsed architectures)")

# thickness: n * dx  (node i occupies [x_i, x_i+dx) in the MATLAB mesh)
def find_block(lines, start_pat, from_idx=0):
    for i in range(from_idx, len(lines)):
        if lines[i].startswith(start_pat):
            return i
    raise KeyError(start_pat)

i_back = find_block(T, "back contact")
i_layer1 = find_block(T, "layer")
i_front = find_block(T, "front contact")
i_work = find_block(T, "working point")
header = T[:i_back]
# the template has three layers; take the last one (poly-Se block with the SRH defect) as the layer template
idx_layers = [i for i, l in enumerate(T) if l == "layer"]
i_lay3 = idx_layers[-1]
layer_tmpl = T[i_lay3:i_front]
back_tmpl = T[i_back:i_layer1]
front_tmpl = T[i_front:i_work]
work_tmpl = T[i_work:]
interface_tmpl = ["interface properties", "", "interfacename : {name}", "Intraband tunneling : 0 (0: no tunneling, 1: with tunneling)",
                  "Relative electron mass :  1.00e+00", "Relative hole mass :  1.00e+00", "", ""]

num_re = re.compile(r"[-+]?\d+\.?\d*(?:[eE][-+]?\d+)?")
def set_vals(lines, key, value, fmt="%.6e"):
    """replace first, 6th and 7th numeric entries of the line starting with `key :`"""
    for i, l in enumerate(lines):
        if l.startswith(key + " :") or l.startswith(key + ":"):
            head, tail = l.split(":", 1)
            toks = list(num_re.finditer(tail))
            new = tail
            for j in sorted([0, 5, 6], reverse=True):
                if j < len(toks):
                    m = toks[j]
                    new = new[:m.start()] + (fmt % value) + new[m.end():]
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
    set_vals(L, "eps", r["eps"], "%.6f")
    set_vals(L, "chi", r["chi"], "%.6f")
    set_vals(L, "Eg", r["Eg"], "%.6f")
    set_vals(L, "Nc", r["Nc"] * 1e6)
    set_vals(L, "Nv", r["Nv"] * 1e6)
    set_vals(L, "mu_n", r["muN"] * 1e-4)
    set_vals(L, "mu_p", r["muP"] * 1e-4)
    set_vals(L, "Na(uniform)", r["Na"] * 1e6)
    set_vals(L, "Nd(uniform)", r["Nd"] * 1e6)
    # SRH: tau = 1/(sigma*vth*Nt); sigma = 1e-19 m2 (1e-15 cm2), vth = 1e5 m/s (template value) -> Nt = 1/(tau*sigma*vth)
    sigma = 1e-19; vth = 1e5; Nt = 1.0 / (r["tau"] * sigma * vth)
    set_vals(L, "Nt(uniform)", Nt)
    set_scalar(L, "Nleft", Nt); set_scalar(L, "Nright", Nt)
    for i, l in enumerate(L):
        if l.startswith("sigma_n"): L[i] = "sigma_n : 1.000e-19\t[m^2]"
        if l.startswith("sigma_p"): L[i] = "sigma_p : 1.000e-19\t[m^2]"
        if l.startswith("Reference for defect energy"):
            L[i] = "Reference for defect energy :  1 ; 0 = above EV (SCAPS<2.7); 1 = above Ei; 2 = below EC"
        if l.startswith("Et :"):
            L[i] = "Et :   0.000\t[eV]"
        if l.startswith("Ekar"):
            L[i] = "Ekar :  0.0000 [eV]"
    return L

out = []
out += ["> SCAPS 3.3.12 cross-check of the paper architecture %s (exported from the MATLAB model)." % tag,
        "> SRH only (midgap neutral defect, tau = 1/(sigma vth Nt)), ideal heterointerfaces, generation from file (MATLAB G(x)).",
        "> Generated by external_validation/scaps_build.py - do not edit.",
        "convergence"] if False else []
hdr = list(header)
# replace the free-text comment lines (those starting with '>') with our own
hdr = [l for l in hdr if not l.startswith(">")]
hdr = ["> SCAPS 3.3.12 cross-check of the paper architecture %s (exported from the MATLAB model)." % tag,
       "> SRH only (midgap neutral defect), ideal heterointerfaces, generation from file (MATLAB G(x))."] + hdr
out += hdr
# left contact ("back contact" block of the file = left side): ideal for both carriers

kT_eV = 0.025852
def flat_phi(r):
    """work function (eV) that makes the contact flat-band for the doping of layer record r (SCAPS ignores the flatband flag in scripts)"""
    Nd_, Na_ = r["Nd"], r["Na"]
    if Nd_ >= Na_:
        return r["chi"] + kT_eV * np.log(r["Nc"] / max(Nd_ - Na_, 1.0))
    return r["chi"] + r["Eg"] - kT_eV * np.log(r["Nv"] / (Na_ - Nd_))
bl = list(back_tmpl)
set_scalar(bl, "Sn", 1e5); set_scalar(bl, "Sp", 1e5); set_scalar(bl, "Fi_m", flat_phi(layer_info[0][3]), "%.6f")
for i, l in enumerate(bl):
    if l.startswith("flatband"): bl[i] = "flatband : 1 (1: flatband contact imposed; 0: fixed Phi_m value)"
out += bl
names = {1: "L1", 2: "L2", 3: "L3", 4: "L4", 5: "L5", 6: "L6"}
for n, (k, cnt, dx, r) in enumerate(layer_info):
    d_m = cnt * dx
    out += make_layer(names[int(k)], d_m, r)
    if n < len(layer_info) - 1:
        out += [s.format(name=f"{names[int(k)]} / {names[int(layer_info[n+1][0])]}") for s in interface_tmpl]
fl = list(front_tmpl)
set_scalar(fl, "Sn", meta["Sback"] * 1e-2); set_scalar(fl, "Sp", 1e5); set_scalar(fl, "Fi_m", flat_phi(layer_info[-1][3]), "%.6f")
for i, l in enumerate(fl):
    if l.startswith("flatband"): fl[i] = "flatband : 1 (1: flatband contact imposed; 0: fixed Phi_m value)"
out += fl
out += work_tmpl
open(os.path.join(work, f"{tag}.def"), "w", encoding="latin-1").write("\n".join(out))

# generation file: x [m], g [m^-3 s^-1]; node centres
xc = []
for k, cnt, dx, r in layer_info:
    pass
dxn = np.empty_like(x)
for k, cnt, dx, r in layer_info:
    sel = nodes["layer"] == k
    dxn[sel] = dx if dx is not None else 0
xcen = x + 0.5 * dxn
g = nodes["G"] * 1e6     # cm^-3 s^-1 -> m^-3 s^-1
Ltot = sum(cnt * dx for k, cnt, dx, r in layer_info)
# SCAPS needs generation-file points exactly on the layer interfaces (otherwise the Newton solver fails at V = 0 in
# a seemingly random way depending on the number of points); a step is written as two points EPS apart.
EPS = 1e-3   # um
pts = [(0.0, g[0])]
pts_layers = [np.where(nodes["layer"] == k)[0] for k in layers]
for n_, idx in enumerate(pts_layers):
    pts += [(xcen[i] * 1e6, g[i]) for i in idx]
    if n_ < len(pts_layers) - 1:
        bnd = (x[pts_layers[n_ + 1][0]]) * 1e6          # first node of the next layer = interface
        pts += [(bnd - EPS, g[idx[-1]]), (bnd + EPS, g[pts_layers[n_ + 1][0]])]
pts.append((Ltot * 1e6, g[-1]))
pts.sort(key=lambda p: p[0])
with open(os.path.join(work, f"{tag}.gen"), "w", newline="") as f:   # newline="" keeps the explicit CRLF (no CR CR LF)
    f.write("Generation profile exported from MATLAB\r\n\r\nx (um)\tGeh (#/m3.s)\r\n")
    for a, b in pts:
        f.write("%.8f\t%.5E\r\n" % (a, b))

scr = f"""// SCAPS cross-check {tag}
set errorhandling.overwritefile
set script_display_mode.fully_suppressed
set batch_display_mode.suppressed
set quitscript.quitSCAPS
load definitionfile {os.path.join(work, tag + '.def')}
clear all
clear actions
action iv.doiv
action iv.startV 0.0
action iv.stopV 1.60
action iv.increment 0.02
action iv.stopaftervoc 0
load generationfile {os.path.join(work, tag + '.gen')}
action generationfromfile
action light
calculate singleshot
save results.iv {os.path.join(work, tag + '_scaps.iv')}
"""
open(os.path.join(work, f"{tag}.script"), "w").write(scr)
print("built", tag, "layers:", [(int(k), cnt, dx) for k, cnt, dx, r in layer_info], "L=%.4e m" % Ltot)
