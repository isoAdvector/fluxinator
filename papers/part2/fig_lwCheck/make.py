# Copyright (C) 2026 Johan Roenby
# SPDX-License-Identifier: GPL-3.0-or-later
"""Is the Lax-Wendroff correction worth keeping? fluxinator (explicit + LW + interfaceUpwind,
run/progressiveWave/finalSolver_maxCo) against the same cases without LW (explicit +
interfaceUpwind, run/progressiveWave/finalSolverNoLW), 190x40 and 380x80 at maximum Courant
numbers 0.5, 0.2 and 0.1, 30 s, with the metrics of finalSolver.py; and the fixed time step 0.02 s
on 190x40 (convection/explicitLWIU_190x40_dt0.02 with LW, finalSolverNoLW/f1_190x40_dt0.02 without).
Part II, Figure fig:lwCheck. Output: lwCheck.pdf, lwCheck.txt. Usage: python3 make.py"""
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

import os
import numpy as np
import finalSolver as F
from common import G, WAVE, CONV, write, plt, save, FULL

FIN = WAVE
CASES = {"0014": (190, 0.5), "0017": (190, 0.2), "0020": (190, 0.1),
         "0026": (380, 0.5), "0029": (380, 0.2), "0032": (380, 0.1)}
VARIANTS = [("LW + interfaceUpwind", "finalSolver_maxCo"), ("interfaceUpwind", "finalSolverNoLW")]
KEYS = [("dH", "{:+6.1f}"), ("dc", "{:+6.2f}"), ("Umax", "{:5.2f}"), ("Uair", "{:5.2f}"),
        ("emax", "{:5.1f}"), ("r0", "{:5.2f}"), ("ulow", "{:+6.3f}"), ("steps", "{:6d}")]


def last_time(d):
    f = os.path.join(d, "postProcessing", "interfaceHeight1", "0", "height.dat")
    if not os.path.exists(f):
        return 0.0
    return float(open(f, "rb").read()[-300:].decode(errors="ignore").strip().splitlines()[-1].split()[0])


RES = {}
L = ["LW check: explicit convection + interfaceUpwind with and without the Lax-Wendroff correction.",
     "Progressive wave of Part I, 30 s, Courant-limited; metrics as in finalSolver.txt (script finalSolver.py):",
     "dH wave height error after 15 periods [%], dc mean celerity error [%], Umax median domain max|U|/u_crest,",
     "Uair median max|U| in air/u_crest, emax crest-column speed error 28-30 s [% of crest-surface speed],",
     "r0 u/u_exact top water cell, ulow mean current at depth/u_crest, steps time steps.",
     "Runs: run/progressiveWave/finalSolver_maxCo/<case> (with LW), finalSolverNoLW/<case> (without).", "",
     f"{'mesh':8s} {'maxCo':>5s} {'convection':22s} " + " ".join(f"{k:>6s}" for k, _ in KEYS)]
for case, (nx, co) in CASES.items():
    for v, sub in VARIANTS:
        d = os.path.join(FIN, sub, case)
        if not F.finished(d):
            L.append(f"{nx:4d}x{nx*4//19:<3d} {co:5g} {v:22s} not finished (t = {last_time(d):.1f} s)")
            continue
        r = F.analyse(("maxCo", case, v, nx, nx*4//19, co, d))
        r["emax"] *= 100
        RES[(nx, co, v)] = r
        L.append(f"{nx:4d}x{nx*4//19:<3d} {co:5g} {v:22s} " +
                 " ".join(fmt.format(r[k]).rjust(6) for k, fmt in KEYS))

L += ["", "Fixed time step 0.02 s, 190x40 (max Courant number 0.65-0.83): time reached of 10 s"]
for v, d in [("LW + interfaceUpwind", os.path.join(CONV, "explicitLWIU_190x40_dt0.02")),
             ("interfaceUpwind", os.path.join(FIN, "finalSolverNoLW", "f1_190x40_dt0.02"))]:
    # the older runs have log.interIsoGoam or log.interIsoGoamExt
    lg = next((p for p in (os.path.join(d, "log." + a) for a in ("fluxinator", "interIsoGoam", "interIsoGoamExt"))
               if os.path.exists(p)), os.path.join(d, "log.fluxinator"))
    t = [float(l.split()[2]) for l in open(lg, errors="ignore") if l.startswith("Time = ")] if os.path.exists(lg) else [0]
    end = "\nEnd\n" in open(lg, errors="ignore").read() if os.path.exists(lg) else False
    L.append(f"  {v:22s} last time {t[-1]:6.2f} s  {'finished' if end else 'stopped (diverged or running)'}")
write("lwCheck.txt", L)
print("\n".join(L))

# figure: wave height error and crest-column speed error against the maximum Courant number
fig, axes = plt.subplots(1, 2, figsize=(FULL, 2.7))
COLM = {190: "#2a78d6", 380: "#eb6834"}
for ax, k, lab in [(axes[0], "dH", "(a) Wave height error after 15 periods [%]"),
                   (axes[1], "emax", "(b) Crest-column speed error [%]")]:
    for nx in (190, 380):
        for v, ls, mk in [("LW + interfaceUpwind", "-", "o"), ("interfaceUpwind", "--", "s")]:
            pts = sorted((co, RES[(n, co, vv)][k]) for (n, co, vv) in RES if n == nx and vv == v)
            if pts:
                x, y = zip(*pts)
                ax.plot(x, y, ls, marker=mk, color=COLM[nx], ms=4, mfc="white" if ls == "--" else COLM[nx],
                        label=f"{nx}x{nx*4//19}, {'with' if v.startswith('LW') else 'without'} LW")
    ax.set_xscale("log"); ax.set_xticks([0.1, 0.2, 0.5]); ax.set_xticklabels(["0.1", "0.2", "0.5"])
    ax.minorticks_off(); ax.set_xlabel("max Courant number"); ax.set_title(lab, fontsize=9, loc="left")
    ax.axhline(0, color="#1f1f1e", lw=0.6)
axes[1].legend(fontsize=8, loc="upper left")
fig.tight_layout()
save(fig, "lwCheck")
