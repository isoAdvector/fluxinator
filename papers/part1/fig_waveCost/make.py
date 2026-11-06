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

# Cost under a Courant limit
MESHCOL = {"95": "#c4c4c0", "190": "#909090", "380": "#5c5c5a", "760": "#1f1f1e"}
by = {(r["application"], r["nx"], float(r["maxCo"])): r for r in limit}
fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.7), sharey=True)
for ax, key, title in [(axes[0], "steps", "(a) number of time steps"),
                       (axes[1], "clockTime", "(b) wall-clock time")]:
    for nx, label in MESHES:
        cos = sorted({c for (a, n, c) in by if n == nx}, reverse=True)
        for other, ls, mfc in [("interFoam", "-", None), ("interIsoFoam", "--", "white")]:
            y = [100*(float(by[("interIsoGoam", nx, c)][key])/float(by[(other, nx, c)][key]) - 1)
                 for c in cos]
            ax.plot(cos, y, ls, color=MESHCOL[nx], marker="o", mfc=mfc or MESHCOL[nx],
                    label=label if other == "interFoam" else None)
    ax.axhline(0, color="#1f1f1e", lw=0.6)
    ax.set_xscale("log")
    ax.set_xticks([0.05, 0.1, 0.2, 0.5])
    ax.set_xticklabels(["0.05", "0.1", "0.2", "0.5"])
    ax.minorticks_off()
    ax.grid(True, lw=0.4, alpha=0.5)
    ax.set_xlabel("maximum Courant number (limit)")
    ax.set_title(title, fontsize=9, loc="left")
axes[0].set_ylabel("change, interIsoGoam vs other [%]")
h = [plt.Line2D([], [], color=MESHCOL[n], marker="o", label=l) for n, l in MESHES]
h += [plt.Line2D([], [], color="#6b6b6b", ls="-", marker="o", label="vs interFoam"),
      plt.Line2D([], [], color="#6b6b6b", ls="--", marker="o", mfc="white", label="vs interIsoFoam")]
fig.legend(handles=h, loc="upper center", ncol=6, bbox_to_anchor=(0.5, 1.0))
fig.tight_layout(rect=(0, 0, 1, 0.9))
fig.savefig(os.path.join(HERE, "waveCost.pdf"))
fig.savefig(os.path.join(HERE, "waveCost.png"), dpi=150)
