# Copyright (C) 2026 Johan Roenby
# SPDX-License-Identifier: GPL-3.0-or-later
"""Part I, Figure fig:waveCost: runs under a Courant limit, change of interIsoGoam relative to interFoam
and interIsoFoam in the number of time steps and the wall-clock time (see ../waveRuns.py).
Output: waveCost.pdf. Usage: python3 make.py"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
from waveRuns import *  # noqa: E402,F401,F403
HERE = os.path.dirname(os.path.abspath(__file__))

# Cost under a Courant limit. Colour: the solver compared with (as in the other wave figures);
# marker and shade: the mesh, light to dark from coarse to fine.
import matplotlib.colors as mc  # noqa: E402
SHADE = {"95": 0.45, "190": 0.65, "380": 0.85, "760": 1.0}
MARK = {"95": "^", "190": "s", "380": "o", "760": "D"}


def col(other, nx):
    c = mc.to_rgb(COLOUR[other])
    s = SHADE[nx]
    return tuple(1 - s*(1 - x) for x in c) if s < 1 else tuple(0.7*x for x in c)


by = {(r["application"], r["nx"], float(r["maxCo"])): r for r in limit}
fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.1), sharey=True)
for ax, key, title in [(axes[0], "steps", "(a) number of time steps"),
                       (axes[1], "clockTime", "(b) wall-clock time")]:
    for nx, label in MESHES:
        cos = sorted({c for (a, n, c) in by if n == nx}, reverse=True)
        for other, ls, mfc in [("interFoam", "-", None), ("interIsoFoam", "--", "white")]:
            y = [100*(float(by[("interIsoGoam", nx, c)][key])/float(by[(other, nx, c)][key]) - 1)
                 for c in cos]
            ax.plot(cos, y, ls, color=col(other, nx), marker=MARK[nx], ms=4.5, mfc=mfc or col(other, nx),
                    label=label if other == "interFoam" else None)
    ax.axhline(0, color="#1f1f1e", lw=0.6)
    ax.set_xscale("log")
    ax.set_xticks([0.05, 0.1, 0.2, 0.5])
    ax.set_xticklabels(["0.05", "0.1", "0.2", "0.5"])
    ax.minorticks_off()
    ax.grid(True, lw=0.4, alpha=0.5)
    ax.set_xlabel("maximum Courant number (limit)")
    ax.set_title(title, loc="left")
axes[0].set_ylabel("change, interIsoGoam vs other [%]")
h = [plt.Line2D([], [], color="#5c5c5a", ls="none", marker=MARK[n], label=l) for n, l in MESHES]
h += [plt.Line2D([], [], color=COLOUR["interFoam"], ls="-", label="vs interFoam"),
      plt.Line2D([], [], color=COLOUR["interIsoFoam"], ls="--", label="vs interIsoFoam (open)")]
fig.legend(handles=h, loc="upper center", ncol=3, bbox_to_anchor=(0.5, 1.0))
fig.tight_layout(rect=(0, 0, 1, 0.83))
fig.savefig(os.path.join(HERE, "waveCost.pdf"))
fig.savefig(os.path.join(HERE, "waveCost.png"), dpi=150)
