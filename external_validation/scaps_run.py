"""Run one SCAPS-1D 3.3.12 script from Python and report the error log.
Usage (in code): from scaps_run import run; run('name', ['action dark', ...], defname='RBSF-GaP', vstop=0.3)
Writes <work>/<name>.script, runs SCAPS (timeout, kills only its own PID), returns (exit code, error-log text, iv path or None).
"""
import os, subprocess, time

# Portable paths (2026-11-07): the work folder and the error log follow the current user unless SCAPS_WORK is set.
WORK = os.environ.get("SCAPS_WORK", os.path.join(os.path.expanduser("~"), "scaps_xcheck"))
EXE = os.environ.get("SCAPS_EXE", r"C:\Program Files (x86)\Scaps3312\scaps3312.exe")   # SCAPS-1D 3.3.12 executable (not distributed)
LOG = os.path.join(os.environ.get("LOCALAPPDATA", ""), "VirtualStore", "Program Files (x86)", "Scaps3312", "SCAPSErrorLogFile.log")


def run(name, body, defname="RBSF-GaP", vstart=0.0, vstop=0.3, dv=0.1, stopvoc=1, timeout=300, post=()):
    iv = os.path.join(WORK, name + ".iv")
    lines = ["set errorhandling.overwritefile", "set script_display_mode.fully_suppressed",
             "set batch_display_mode.suppressed", "set quitscript.quitSCAPS",
             "load definitionfile " + os.path.join(WORK, defname + ".def"),
             "clear all", "clear actions", "action iv.doiv",
             "action iv.startV %g" % vstart, "action iv.stopV %g" % vstop, "action iv.increment %g" % dv,
             ("action iv.stopaftervoc" if stopvoc else "action iv.continueafterVoc")] + list(body) + \
            ["calculate singleshot", "save results.iv " + iv] + list(post)
    sp = os.path.join(WORK, name + ".script")
    open(sp, "w").write("\n".join(lines) + "\n")
    if os.path.exists(iv): os.remove(iv)
    if os.path.exists(LOG): os.remove(LOG)
    p = subprocess.Popen([EXE, sp], cwd=os.path.dirname(EXE))
    t0 = time.time()
    try:
        rc = p.wait(timeout=timeout)
    except subprocess.TimeoutExpired:
        p.kill(); rc = "TIMEOUT"
    log = open(LOG, encoding="latin-1").read().replace("\n", " | ") if os.path.exists(LOG) else ""
    return rc, log, (iv if os.path.exists(iv) else None), time.time() - t0


if __name__ == "__main__":
    gen = os.path.join(WORK, "RBSF-GaP.gen")
    for nm, body in [("tD", ["action dark"]), ("tL", ["action light"]),
                     ("tG", ["action light", "load generationfile " + gen, "action generation_specified"])]:
        print(nm, run(nm, body))
