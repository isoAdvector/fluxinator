# Copyright (C) 2026 Johan Roenby
# SPDX-License-Identifier: GPL-3.0-or-later
"""Part I, Figure fig:waveCost: runs under a Courant limit, change of interIsoGoam relative to interFoam
and interIsoFoam in the number of time steps and the wall-clock time (see ../waveRuns.py).
The wall-clock time of a run is its number of steps times the time per step measured in the
timing runs (run/progressiveWave/timing_partI, data/timing_partI.csv). The change of
interIsoGoam relative to another solver is formed in each repeat in which both ran at the
same time, and the median over these repeats is shown (two on the coarser meshes, three on
760x160).
Output: waveCost.pdf and waveCost.txt (the numbers, per repeat, and the cost of the second
layer). Usage: python3 make.py"""
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

# Time per step [ms] per (application, nx, maxCo, nLayers) and repeat
import csv  # noqa: E402
from statistics import median  # noqa: E402
ms = {}
for r in csv.DictReader(open(os.path.join(PD, "timing_partI.csv"))):
    ms.setdefault((r["application"], r["nx"], float(r["maxCo"]), r["nLayers"]), {})[r["repeat"]] = \
        float(r["msPerStep"])
def wall_changes(nx, co, other, a="interIsoGoam", nL=("2", "2")):
    """Change [%] of the wall-clock time of a relative to other, per shared repeat"""
    ga, go = ms[(a, nx, co, nL[0])], ms[(other, nx, co, nL[1])]
    steps = int(by[(a, nx, co)]["steps"])/int(by[(other, nx, co)]["steps"]) if other != a else 1
    return [100*(steps*ga[k]/go[k] - 1) for k in sorted(set(ga) & set(go))]


txt = ["Runs under a Courant limit: change of interIsoGoam relative to interFoam (F) and",
       "interIsoFoam (IF) [%] in the number of steps and the wall-clock time (steps of the full",
       "run times the median time per step of the timing runs; each repeat alone in brackets),",
       "and the median time per step [ms].",
       f"{'mesh':>5s} {'maxCo':>6s} {'steps F':>8s} {'steps IF':>9s} {'wall F':>28s} {'wall IF':>28s}"
       f" {'ms F':>8s} {'ms IF':>8s} {'ms IG':>8s}"]
for nx, _ in MESHES:
    for co in sorted({c for (a, n, c) in by if n == nx}, reverse=True):
        g = by[("interIsoGoam", nx, co)]
        line = f"{nx:>5s} {co:6.2f}"
        for o in SOLVERS[:2]:
            line += f" {100*(int(g['steps'])/int(by[(o, nx, co)]['steps']) - 1):8.1f}"
        for o in SOLVERS[:2]:
            per = wall_changes(nx, co, o)
            line += f" {median(per):6.1f} [" + ",".join(f"{x:5.1f}" for x in per) + "]" + " "*6*(3 - len(per))
        line += "".join(f" {median(ms[(a, nx, co, '2')].values()):8.1f}" for a in SOLVERS)
        txt.append(line)
txt += ["", "Second layer: time per step of interIsoGoam with nLayers 2 relative to nLayers 1 in the",
        "same round, 380x80 [%], per repeat"]
for co in (0.5, 0.2, 0.1, 0.05):
    per = wall_changes("380", co, "interIsoGoam", nL=("2", "1"))
    txt.append(f"  maxCo {co:4.2f}: " + ", ".join(f"{x:5.1f}" for x in per))
open(os.path.join(HERE, "waveCost.txt"), "w").write("\n".join(txt) + "\n")
fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.1), sharey=True)
for ax, key, title in [(axes[0], "steps", "(a) number of time steps"),
                       (axes[1], "clockTime", "(b) wall-clock time")]:
    for nx, label in MESHES:
        cos = sorted({c for (a, n, c) in by if n == nx}, reverse=True)
        for other, ls, mfc in [("interFoam", "-", None), ("interIsoFoam", "--", "white")]:
            y = [median(wall_changes(nx, c, other)) if key == "clockTime" else
                 100*(float(by[("interIsoGoam", nx, c)][key])/float(by[(other, nx, c)][key]) - 1)
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
