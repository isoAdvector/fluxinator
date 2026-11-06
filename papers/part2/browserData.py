"""Data of the selection figures (selection.py) and the results browser: the Part 1 solvers rerun
with ddtCorr off and fluxinator (run/progressiveWave/finalSolver_*), progressive wave of
Part 1, 30 s.
Same case numbers in all sets; only finished runs (height.dat reaches 30 s).

Per case (both scans): dH(t) per wave (analyseScan.load) and dHfit, celerity error dcMean,
domain max|U|(t)/u_crest (fieldMinMax, binned maxima), steps (rows of fieldMinMax.dat), clock
time (ExecutionTime of the solver log; missing for the three Part 1 runs whose log was lost).
Courant-limited scan only (fields every 0.2 s): max|U| in the air (alpha < 0.01) per snapshot,
emax, r0, ulow over 28-30 s as in finalSolver.py, and at t = 10, 20, 30 s the crest-column
profile u_x(z) with the exact profile at the solver's crest phase (kinematics.snapshot) and the
surface elevation h(x) relative to the crest position with the exact elevation.
Writes data/browser.json. Usage: python3 browserData.py [-j N]"""
import glob
import json
import os
import re
import sys
from multiprocessing import Pool
import numpy as np
from common import G, K, WAVE, DATA
import finalSolver as FS

FLUME = WAVE
sys.path.insert(0, FLUME)
from analyseScan import load, metrics, U_CREST  # noqa: E402

FINAL = FS.FINAL
SNAPS = ("10", "20", "30")


def cases():
    return FS.cases()


def r4(a):
    return [float(f"{v:.4g}") for v in np.asarray(a, dtype=float)]


def clock(d):
    for f in sorted(glob.glob(os.path.join(d, "log.inter*")) + glob.glob(os.path.join(d, "log.fluxinator"))):
        ex = re.findall(r"^ExecutionTime = (\S+)", open(f, errors="ignore").read(), re.M)
        if ex:
            return float(ex[-1])
    return None


def binned_max(t, u, n=300):
    if len(t) <= n:
        return t, u
    edges = np.linspace(t[0], t[-1], n + 1)
    idx = np.clip(np.searchsorted(edges, t, side="right") - 1, 0, n - 1)
    tb, ub = [], []
    for i in range(n):
        m = idx == i
        if m.any():
            j = np.argmax(u[m])
            tb.append(t[m][j]); ub.append(u[m][j])
    return np.array(tb), np.array(ub)


def surface(d, t, M, W):
    A = K.foamio.read_vol(f"{d}/{t}/alpha.water").reshape(M.nz, M.nx)
    h = A.sum(0)*M.dz                       # elevation above the bed
    ic = int(np.argmax(h)); sc, _ = K._parab(h + M.z1, ic, 1)
    xr = (np.arange(M.nx) + 0.5 - sc)*M.dx
    xr -= M.L*np.round(xr/M.L)
    o = np.argsort(xr)
    return xr[o]/M.L, (h[o] - W.d)/W.H


def fields(d):
    M, W = K.Mesh(d), K.Wave(d)
    ts = [t for t in K.times(d) if float(t) > 2.0 - 1e-9]
    tair, air, late = [], [], []
    for t in ts:
        A = K.foamio.read_vol(f"{d}/{t}/alpha.water")
        U = K.foamio.read_vol(f"{d}/{t}/U").reshape(-1, 3)
        tair.append(float(t)); air.append(np.sqrt((U[A < 0.01]**2).sum(1)).max()/U_CREST)
        if float(t) > 28.0 - 1e-9:
            late.append(K.snapshot(d, t, M, W, None))
    eta = W.eta(0.0)[0]
    low = []
    for r in late:
        zz, un, ue = r["prof"]
        m = zz - M.z1 < 0.5*eta
        low.append(np.mean(un[m] - ue[m])/U_CREST)
    snaps = {}
    for t in SNAPS:
        if not os.path.isdir(f"{d}/{t}"):
            continue
        s = K.snapshot(d, t, M, W, None)
        z, un, ue = s["prof"]
        x, h = surface(d, t, M, W)
        snaps[t] = dict(z=r4((z - M.z1)/eta), u=r4(un/U_CREST), ue=r4(ue/U_CREST),
                        x=r4(x), h=r4(h))
    return dict(tair=r4(tair), Uair=r4(air), UairMed=float(np.median(air)),
                emax=float(100*np.mean([r["emax"] for r in late])),
                r0=float(np.mean([r["r0"] for r in late])), ulow=float(np.mean(low)),
                snaps=snaps)


REF = None


def analyse(c):
    global REF
    scan, case, s, nx, nz, par, d = c
    if REF is None:
        REF = load(os.path.join(WAVE, "solvers_fixedDt"), "0036")
    r = load(os.path.dirname(d), case)
    m = metrics(r, REF, REF["t"][-1])
    Hx = REF["H"][:, 1].mean()
    tu, umax = np.loadtxt(os.path.join(d, "postProcessing/fieldMinMax1/0/fieldMinMax.dat"),
                          comments="#", unpack=True, usecols=(0, 2))
    tb, ub = binned_max(tu, umax/U_CREST)
    row = dict(scan=scan, case=case, solver=s, nx=nx, nz=nz, par=par,
               dHfit=m.get("dHfit"), dc=m.get("dcMean"), steps=int(len(tu)), clock=clock(d),
               UmaxMed=float(np.median(umax[tu > 2])/U_CREST),
               tH=r4(r["H"][:, 0]), dH=r4((r["H"][:, 1] - Hx)/Hx*100),
               tU=r4(tb), Umax=r4(ub))
    if scan == "maxCo":
        row.update(fields(d))
    return row


def exact():
    """Exact surface and crest profile, from any fluxinator case (same wave everywhere)."""
    d = os.path.join(WAVE, "finalSolver_maxCo", "0017")
    W = K.Wave(d)
    th = np.linspace(-np.pi, np.pi, 401)
    eta = W.eta(0.0)[0]
    z = np.linspace(0, eta, 101)
    return dict(x=r4(th/(2*np.pi)), h=r4((W.eta(th) - W.d)/W.H), z=r4(z/eta),
                u=r4(W.uv(0.0*z, z)[0]/U_CREST), H=W.H, d=W.d, L=W.L, T=2.0,
                uCrest=U_CREST)


def main(nproc=8):
    C = cases()
    print(f"{len(C)} finished runs", flush=True)
    with Pool(nproc) as p:
        rows = p.map(analyse, C, chunksize=1)
    out = dict(exact=exact(), cases=rows)
    with open(os.path.join(DATA, "browser.json"), "w") as f:
        json.dump(out, f, separators=(",", ":"))
    print("written", os.path.getsize(os.path.join(DATA, "browser.json"))//1024, "kB")


if __name__ == "__main__":
    main(int(sys.argv[sys.argv.index("-j") + 1]) if "-j" in sys.argv else 8)
