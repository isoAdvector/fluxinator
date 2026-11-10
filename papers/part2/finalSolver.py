# Copyright (C) 2026 Johan Roenby
# SPDX-License-Identifier: GPL-3.0-or-later
"""fluxinator (explicit convection, Gauss laxWendroff ... interfaceUpwind, ddtCorr off;
run/progressiveWave/finalSolver_*) against the Part 1 solvers rerun with ddtCorr off
(run/progressiveWave/solvers_*_ddtCorrOff), on the progressive wave of Part 1 (H 0.125 m, T 2 s,
d 0.4 m, one wavelength, periodic, 30 s = 15 periods). Same case numbers in both sets.
Only finished runs are used (Part 1 reruns of the 760x160 mesh are partly still running).

Per case:
  dH      wave height error after 15 periods, dHfit of run/progressiveWave/analyseScan.py [%]
  dc      mean celerity error over the run, dcMean of analyseScan.py [%]
  Umax    median over t > 2 s of the domain max|U| (fieldMinMax, every step) / u_crest
  steps   number of time steps (rows of fieldMinMax.dat)
Courant-limited scan only (fields written every 0.2 s):
  Uair    median over the snapshots t > 2 s of max|U| in the air (alpha.water < 0.01) / u_crest
  emax, r0  crest column against the exact input wave, mean over 28-30 s (kinematics.py):
          emax = max over wet cells of | |U| - |U_exact| | / exact crest-surface speed,
          r0 = u_x/u_x,exact in the top full water cell,
          ulow = mean of (u_x - u_x,exact)/u_crest over the crest-column cells below half the
          exact crest elevation (a uniform offset: mean current at depth)
  profile crest column u_x(z) at t = 30 s and the exact profile at the same phase
Writes fig_finalOverview/finalSolver.txt, finalSolver.csv and finalOverview.pdf (run by its make.py).
Also imported by fig_lwCheck/make.py and browserData.py."""
import csv
import os
import sys
from multiprocessing import Pool
import numpy as np
from common import G, K, WAVE, plt, save, write, outpath, FULL, DATA, STYLE

FLUME = WAVE
sys.path.insert(0, WAVE)
from analyseScan import load, metrics, U_CREST  # noqa: E402

# plotCrestHeight (imported by analyseScan) sets its own style
plt.rcdefaults()
plt.rcParams.update(STYLE)

SCANS = {"fixedDt": ("solvers_fixedDt_ddtCorrOff", "finalSolver_fixedDt", "solvers_fixedDt.csv", "deltaT"),
         "maxCo": ("solvers_maxCo_ddtCorrOff", "finalSolver_maxCo", "solvers_maxCo.csv", "maxCo")}
FINAL = "explicit + LW + interfaceUpwind"
SOLVERS = ["interFoam", "interIsoFoam", "interIsoGoam", FINAL]
COLOUR = {"interFoam": "#2a78d6", "interIsoFoam": "#eb6834", "interIsoGoam": "#1baf7a",
          FINAL: "#1f1f1e"}
MARKER = {"interFoam": "s", "interIsoFoam": "^", "interIsoGoam": "o", FINAL: "D"}
MESHES = [190, 380, 760]


def cases():
    """(scan, case number, solver, nx, nz, dt or maxCo, case directory) of every finished run."""
    out = []
    for scan, (off, fin, f, par) in SCANS.items():
        for r in csv.DictReader(open(os.path.join(WAVE, "data", f))):
            app = r["application"]
            if app not in SOLVERS:
                continue
            dirs = [(app, os.path.join(WAVE, off, r["case"]))]
            if app == "interIsoGoam":
                dirs.append((FINAL, os.path.join(WAVE, fin, r["case"])))
            for s, d in dirs:
                if finished(d):
                    out.append((scan, r["case"], s, int(r["nx"]), int(r["nz"]), float(r[par]), d))
    return out


def finished(d):
    """height.dat reaches 30 s; logs of three reruns were lost, so the log is not used."""
    f = os.path.join(d, "postProcessing", "surfacePosition1", "0", "height.dat")
    if not os.path.exists(f):
        return False
    last = open(f, "rb").read()[-300:].decode(errors="ignore").strip().splitlines()[-1].split()
    return float(last[0]) > 29.99


def fields(d):
    """Air speed, crest-column errors and the final profile from the written fields."""
    M, W = K.Mesh(d), K.Wave(d)
    air, rows = [], []
    ts = [t for t in K.times(d) if float(t) > 2.0 - 1e-9]
    for t in ts:
        A = K.foamio.read_vol(f"{d}/{t}/alpha.water")
        U = K.foamio.read_vol(f"{d}/{t}/U").reshape(-1, 3)
        air.append(np.sqrt((U[A < 0.01]**2).sum(1)).max()/U_CREST)
        if float(t) > 28.0 - 1e-9:
            rows.append(K.snapshot(d, t, M, W, None))
    z, un, ue = rows[-1]["prof"]
    eta = W.eta(0.0)[0]
    low = []
    for r in rows:
        zz, un, ue = r["prof"]
        m = zz - M.z1 < 0.5*W.eta(0.0)[0]
        low.append(np.mean(un[m] - ue[m])/U_CREST)
    return dict(Uair=float(np.median(air)), emax=float(np.mean([r["emax"] for r in rows])),
                ulow=float(np.mean(low)),
                r0=float(np.mean([r["r0"] for r in rows])),
                prof=((z - M.z1)/(eta), un, ue), dz=M.dz/eta)


REF = None


def analyse(c):
    scan, case, s, nx, nz, par, d = c
    global REF
    if REF is None:
        REF = load(os.path.join(WAVE, "solvers_fixedDt"), "0036")
    m = metrics(load(os.path.dirname(d), case), REF, REF["t"][-1])
    tu, umax = np.loadtxt(os.path.join(d, "postProcessing/fieldMinMax1/0/fieldMinMax.dat"),
                          comments="#", unpack=True, usecols=(0, 2))
    row = dict(scan=scan, case=case, solver=s, nx=nx, nz=nz, par=par, dH=m["dHfit"],
               dc=m["dcMean"], Umax=float(np.median(umax[tu > 2])/U_CREST), steps=len(tu),
               avgCo=m["avgMaxCo"])
    if scan == "maxCo":
        row.update(fields(d))
    return row


def table(rows):
    L = ["fluxinator (Part 2) vs Part 1 solvers with ddtCorr off; progressive wave, 30 s. Script: finalSolver.py",
         "dH: wave height error after 15 periods [%]; dc: mean celerity error [%]; Umax: median domain max|U|/u_crest;",
         "Uair: median max|U| in the air (alpha < 0.01)/u_crest; emax: crest-column speed error / exact crest-surface",
         "speed, 28-30 s [%]; r0: u/u_exact, top full water cell under the crest, 28-30 s; ulow: mean (u - u_exact)/u_crest",
         "in the crest column below half the crest elevation, 28-30 s; steps: time steps.",
         "fluxinator = explicit + LW + interfaceUpwind (Part 2), ddtCorr off. '-': run diverged before 30 s; fields of the fixed-dt Part 1 runs were not written (no Uair, emax, r0)."]
    keys = [("dH", "{:+6.1f}"), ("dc", "{:+6.2f}"), ("Umax", "{:5.2f}"), ("Uair", "{:5.2f}"),
            ("emax", "{:5.1f}"), ("r0", "{:5.2f}"), ("ulow", "{:+6.3f}"), ("steps", "{:6d}")]
    idx = {(r["scan"], r["nx"], r["par"], r["solver"]): r for r in rows}
    for scan, (_, _, _, par) in SCANS.items():
        L += ["", f"{scan} ({par})", f"{'mesh':8s} {par:>9s} {'solver':13s} " +
              " ".join(f"{k:>6s}" for k, _ in keys)]
        for nx in [95] + MESHES:
            for p in sorted({r["par"] for r in rows if r["scan"] == scan and r["nx"] == nx}, reverse=True):
                for s in SOLVERS:
                    r = idx.get((scan, nx, p, s))
                    cells = []
                    for k, fmt in keys:
                        v = None if r is None else r.get(k)
                        if v is not None and k == "emax":
                            v = 100*v
                        cells.append(f"{'-':>6s}" if v is None else fmt.format(v).rjust(6))
                    L.append(f"{nx:4d}x{nx*4//19:<3d} {p:9g} {('fluxinator' if s == FINAL else s):13s} " + " ".join(cells))
    return L


def fig_overview(rows):
    """Courant-limited scan: dH, celerity, Uair, emax, ulow against max Co, one column per mesh."""
    R = [r for r in rows if r["scan"] == "maxCo"]
    spec = [("dH", "wave height\nerror [%]", 1), ("dc", "celerity error [%]", 1),
            ("Uair", "max $|\\mathbf{U}|$ in air\n/ $u_\\mathrm{crest}$", 1),
            ("emax", "crest-column\nspeed error [%]", 100),
            ("ulow", "mean velocity error\nat depth / $u_\\mathrm{crest}$", 1)]
    fig, axs = plt.subplots(len(spec), len(MESHES), figsize=(FULL, 7.4), sharex=True, sharey="row")
    for j, nx in enumerate(MESHES):
        for i, (k, lab, f) in enumerate(spec):
            ax = axs[i, j]
            for s in SOLVERS:
                pts = sorted((r["par"], f*r[k]) for r in R if r["nx"] == nx and r["solver"] == s)
                if pts:
                    x, y = zip(*pts)
                    ax.plot(x, y, "-", color=COLOUR[s], marker=MARKER[s],
                            lw=1.6 if s == FINAL else 1.0)
            ax.axhline(1 if k == "Uair" else 0, color="#6b6b6b", lw=0.6)
            ax.set_xscale("log")
            if j == 0:
                ax.set_ylabel(f"({'abcde'[i]}) " + lab)
        axs[0, j].set_title(f"{nx}x{nx*4//19}")
        axs[-1, j].set_xlabel("max Courant number")
        axs[-1, j].set_xticks([0.05, 0.1, 0.2, 0.5])
        axs[-1, j].set_xticklabels(["0.05", "0.1", "0.2", "0.5"])
    axs[0, 0].set_yscale("symlog", linthresh=2)
    axs[0, 0].set_yticks([-30, -10, -2, 0, 2, 10, 30])
    axs[0, 0].set_yticklabels(["-30", "-10", "-2", "0", "2", "10", "30"])
    h = [plt.Line2D([], [], color=COLOUR[s], marker=MARKER[s], label="fluxinator" if s == FINAL else s + ", ddtCorr off")
         for s in SOLVERS]
    fig.legend(handles=h, loc="upper center", ncol=2, bbox_to_anchor=(0.5, 1.0))
    fig.tight_layout(rect=(0, 0, 1, 0.955))
    save(fig, "finalOverview")


def main(nproc=8):
    C = cases()
    print(f"{len(C)} finished runs", flush=True)
    with Pool(nproc) as p:
        rows = p.map(analyse, C, chunksize=1)
    keys = ["scan", "case", "solver", "nx", "nz", "par", "dH", "dc", "Umax", "Uair", "emax", "r0",
            "ulow", "steps", "avgCo"]
    with open(outpath("finalSolver.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, keys, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    write("finalSolver.txt", table(rows))
    fig_overview(rows)


if __name__ == "__main__":
    main(int(sys.argv[sys.argv.index("-j") + 1]) if "-j" in sys.argv else 8)
