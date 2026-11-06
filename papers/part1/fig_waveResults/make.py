# Copyright (C) 2026 Johan Roenby
# SPDX-License-Identifier: GPL-3.0-or-later
"""Part I, Figure fig:waveResults: wave height error and median max|U|/u_crest of all runs against the
time-averaged maximum Courant number (see ../waveRuns.py). Output: waveResults.pdf. Usage: python3 make.py"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
from waveRuns import *  # noqa: E402,F401,F403
HERE = os.path.dirname(os.path.abspath(__file__))

fig, axes = plt.subplots(2, 4, figsize=(7.0, 4.4), sharex=True, sharey="row")
for j, (nx, label) in enumerate(MESHES):
    for i, (key, ref) in enumerate([("dHfit", 0.0), ("Umed", 1.0)]):
        ax = axes[i, j]
        for app in SOLVERS:
            x, y = series(fixed, app, nx, key)
            ax.plot(x, y, "-", color=COLOUR[app], marker=MARKER[app])
            x, y = series(limit, app, nx, key)
            ax.plot(x, y, "--", color=COLOUR[app], marker=MARKER[app], mfc="white")
        ax.axhline(ref, color="#1f1f1e", lw=0.6)
        ax.set_xscale("log")
        ax.set_xticks([0.05, 0.1, 0.2, 0.5])
        ax.set_xticklabels(["0.05", "0.1", "0.2", "0.5"])
        ax.minorticks_off()
        ax.grid(True, lw=0.4, alpha=0.5)
        if i == 0:
            ax.set_title(label, fontsize=9)
axes[0, 0].set_ylabel(r"$\Delta H$ after 15 periods [%]")
axes[1, 0].set_ylabel(r"median $\max|\mathbf{U}|/u_{crest}$")
axes[0, 0].set_ylim(-60, 15)
fig.supxlabel("time-averaged maximum Courant number", fontsize=9, y=0.01)
h = [plt.Line2D([], [], color=COLOUR[a], marker=MARKER[a], label=a) for a in SOLVERS]
h += [plt.Line2D([], [], color="#6b6b6b", ls="-", marker="o", label=r"fixed $\Delta t$"),
      plt.Line2D([], [], color="#6b6b6b", ls="--", marker="o", mfc="white", label="Courant limit")]
fig.legend(handles=h, loc="upper center", ncol=5, bbox_to_anchor=(0.5, 1.0))
fig.tight_layout(rect=(0, 0, 1, 0.94))
fig.savefig(os.path.join(HERE, "waveResults.pdf"))
fig.savefig(os.path.join(HERE, "waveResults.png"), dpi=150)

