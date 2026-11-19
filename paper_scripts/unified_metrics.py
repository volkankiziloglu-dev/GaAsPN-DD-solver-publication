"""One extractor for AM1.5G and indoor J-V curves (Python port of src/pvMetrics.m): J_sc interpolated at V = 0, V_oc by linear interpolation of the zero crossing,
maximum power by a parabola through the three points around the grid maximum (step-independent MPP). Applied to all finite points of a curve."""
import numpy as np


def pv_metrics(V, J, Pin_Wm2):
    V = np.asarray(V, float); J = np.asarray(J, float)
    o = np.argsort(V); V = V[o]; J = J[o]
    Jsc = np.nan
    for k in range(len(V) - 1):
        if V[k] <= 0 <= V[k + 1] and V[k + 1] > V[k]:
            t = (0 - V[k]) / (V[k + 1] - V[k]); Jsc = J[k] + t * (J[k + 1] - J[k]); break
        if V[k] == 0: Jsc = J[k]; break
    Voc = np.nan
    for k in range(len(V) - 1):
        if V[k] >= 0 and J[k] >= 0 and J[k + 1] < 0:
            t = (0 - J[k]) / (J[k + 1] - J[k]); Voc = V[k] + t * (V[k + 1] - V[k]); break
    if not (np.isfinite(Jsc) and np.isfinite(Voc) and Voc > 0 and Jsc > 0):
        return dict(Jsc=Jsc, Voc=Voc, FF=np.nan, PCE=np.nan, Pmax=np.nan)
    Pmax = 0.0; kb = -1
    for k in range(len(V)):
        if V[k] < 0 or V[k] > Voc or J[k] < 0: continue
        P = V[k] * J[k]
        if P > Pmax: Pmax = P; kb = k
    if 0 < kb < len(V) - 1:
        V1, V2, V3 = V[kb - 1], V[kb], V[kb + 1]; P1, P2, P3 = V1 * J[kb - 1], Pmax, V3 * J[kb + 1]
        num = (V2 - V1) ** 2 * (P2 - P3) - (V2 - V3) ** 2 * (P2 - P1)
        den = (V2 - V1) * (P2 - P3) - (V2 - V3) * (P2 - P1)
        if abs(den) > 0:
            Vv = V2 - 0.5 * num / den
            if V1 < Vv < V3:
                Pv = (P1 * (Vv - V2) * (Vv - V3) / ((V1 - V2) * (V1 - V3)) + P2 * (Vv - V1) * (Vv - V3) / ((V2 - V1) * (V2 - V3)) +
                      P3 * (Vv - V1) * (Vv - V2) / ((V3 - V1) * (V3 - V2)))
                if Pv > Pmax: Pmax = Pv
    return dict(Jsc=Jsc, Voc=Voc, FF=Pmax / (Voc * Jsc), PCE=Pmax * 1e4 / Pin_Wm2 * 100, Pmax=Pmax)
