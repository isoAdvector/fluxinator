# Copyright (C) 2026 Johan Roenby
# SPDX-License-Identifier: GPL-3.0-or-later
"""Crest kinematics of the three solvers in the stored adaptive runs on the 380x80 mesh
(run/progressiveWave/solvers_maxCo cases 0027 interFoam, 0028 interIsoFoam, 0029 interIsoGoam, max Co 0.2,
30 s). Horizontal velocity under the crest compared with the exact stream function wave of
the input wave (the initial condition of all runs).

Shared by papers/part1/tab_crest and fig_crestProfiles (one make.py each).
Per time window: dH = wave height change, r0 = u/u_exact in the top full water cell under the
crest, emax = max | |U| - |U_exact| | over the wet cells of the crest column [% of the exact
crest-surface speed]; means over the snapshots of the window. The analysis functions are
those of run/common/hmatch.py and run/common/kinematics.py."""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
WAVE = os.path.join(HERE, "..", "..", "run", "progressiveWave")
sys.path.insert(0, os.path.join(HERE, "..", "..", "run", "common"))
import kinematics as K  # noqa: E402
from hmatch import analyse_input  # noqa: E402

CASES = {s: os.path.join(WAVE, "solvers_maxCo", c) for s, c in
         [("interFoam", "0027"), ("interIsoFoam", "0028"), ("interIsoGoam", "0029")]}
COLOUR = {"interFoam": "#DDAA33", "interIsoFoam": "#BB5566", "interIsoGoam": "#004488"}  # Paul Tol high-contrast
MARKER = {"interFoam": "o", "interIsoFoam": "s", "interIsoGoam": "D"}
WINDOWS = [(8, 10), (18, 20), (28, 30)]
TIMES = ["10", "30"]
sys.path.insert(0, HERE)
from common import paper_style  # noqa: E402
paper_style(scale=5.0/7.0)  # the 7 in wide figures are printed at the text width
plt.rcParams.update({"axes.linewidth": 0.7, "lines.linewidth": 1.2, "lines.markersize": 3.5,
                     "legend.frameon": False, "axes.grid": True, "grid.linewidth": 0.4,
                     "grid.alpha": 0.5})

def table(out):
    """Write the crest kinematics table to the file out."""
    lines = ["Crest kinematics, 380x80 mesh, max Co 0.2 (run/progressiveWave/solvers_maxCo 0027-0029)",
             "dH: wave height change [%]; r0: u/u_exact in the top full water cell under the crest;",
             "emax: max ||U|-|U_exact|| in the crest column [% of exact crest-surface speed];",
             "exact = the input stream function wave. Script: papers/part1/tab_crest/make.py", "",
             f"{'solver':14s} " + "  ".join(f"{f'{a}-{b} s: dH   r0    emax':>26s}" for a, b in WINDOWS)]
    for s, c in CASES.items():
        cells = []
        for w in WINDOWS:
            a = analyse_input(c, w)
            cells.append(f"{100*a['dH']:+6.1f} {a['r0']:5.2f} {100*a['emax']:6.1f} %")
        lines.append(f"{s:14s} " + "  ".join(f"{x:>26s}" for x in cells))
    open(out, "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))



def profiles(out):
    """Draw the crest velocity profiles to out.pdf and out.png."""
    fig, axes = plt.subplots(1, len(TIMES), figsize=(7.0, 3.3), sharey=True)
    for ax, t in zip(axes, TIMES):
        for s, c in CASES.items():
            M, W0 = K.Mesh(c), K.Wave(c)
            snap = K.snapshot(c, t, M, W0, None)
            z, un, _ = snap["prof"]
            ax.plot(un, z, "-", marker=MARKER[s], color=COLOUR[s],
                    label=f"{s} ($\\Delta H$ {100*(snap['H']/W0.H - 1):+.0f} %)")
        ze = np.linspace(M.z1, M.z1 + W0.eta(0.0)[0], 200)
        ax.plot(W0.uv(0.0, ze - M.z1)[0], ze, "--", color="#1f1f1e", lw=1.0, label="exact (input wave)")
        ax.set_title(f"$t = {t}$ s")
        ax.set_xlabel("horizontal velocity under the crest [m/s]")
        ax.set_ylim(-0.42, 0.11)
        ax.legend(loc="lower right")
    axes[0].set_ylabel("z above still water level [m]")
    fig.tight_layout()
    fig.savefig(out + ".pdf")
    fig.savefig(out + ".png", dpi=150)
