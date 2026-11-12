# Copyright (C) 2026 Johan Roenby
# SPDX-License-Identifier: GPL-3.0-or-later
"""Part I, Figure fig:ddtCorr: progressive wave with adjusted time step on the 190x40 and
380x80 meshes, ddtCorr on (run/progressiveWave/solvers_maxCo) and off
(run/progressiveWave/solvers_maxCo_ddtCorrOff). dHfit: wave height error after 15 periods [%]
(analyseScan.py); r0: horizontal velocity in the top full water cell under the crest relative
to the exact input wave, mean over 28-30 s (hmatch.analyse_input, as crestKinematics.py).

Figure: (a), (b) wave height error after 15 periods against the maximum Courant number,
all three solvers, ddtCorr on (solid) and off (dashed); (c) interIsoGoam, r0, both meshes
(light: 190x40). Output: ddtCorr.txt (the numbers) and ddtCorr.pdf. Usage: python3 make.py
"""
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
FLUME = os.path.join(HERE, "..", "..", "..", "run", "progressiveWave")
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "run", "common"))
from hmatch import analyse_input  # noqa: E402

on = {r["case"]: r for r in csv.DictReader(open(os.path.join(FLUME, "data", "solvers_maxCo.csv")))}
off = {r["case"]: r for r in csv.DictReader(open(os.path.join(FLUME, "data", "solvers_maxCo_ddtCorrOff.csv")))}
COS = [0.5, 0.2, 0.1, 0.05]
L = ["ddtCorr on / off, adjusted time step, max Co 0.5 / 0.2 / 0.1 / 0.05. Script: papers/part1/fig_ddtCorr/make.py",
     "dHfit [%] after 15 periods (analyseScan.py); r0 = u/u_exact, top full water cell under the crest, 28-30 s", ""]
for nx in (190, 380):
    for app in ("interFoam", "interIsoFoam", "interIsoGoam"):
        cs = {float(r["maxCo"]): c for c, r in on.items() if int(r["nx"]) == nx and r["application"] == app}
        L.append(f"{nx}x{nx*4//19} {app:13s} dHfit " + "  ".join(
            f"{float(on[cs[co]]['dHfit']):+6.1f} / {float(off[cs[co]]['dHfit']):+6.1f}" for co in COS))
        if app == "interIsoGoam":
            r = []
            for co in COS:
                c = cs[co]
                a = analyse_input(os.path.join(FLUME, "solvers_maxCo", c), (28, 30))["r0"]
                b = analyse_input(os.path.join(FLUME, "solvers_maxCo_ddtCorrOff", c), (28, 30))["r0"]
                r.append(f"{a:6.2f} / {b:6.2f}")
            L.append(f"{nx}x{nx*4//19} {app:13s} r0    " + "  ".join(r))
open(os.path.join(HERE, "ddtCorr.txt"), "w").write("\n".join(L) + "\n")
print("\n".join(L))


# ---- figure
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(HERE, ".."))
from common import paper_style  # noqa: E402
paper_style(scale=5.0/7.0)  # 7 in wide, printed at the text width
plt.rcParams.update({"axes.linewidth": 0.7, "lines.linewidth": 1.2, "lines.markersize": 4,
                     "legend.frameon": False})
COLOUR = {"interFoam": "#DDAA33", "interIsoFoam": "#BB5566", "interIsoGoam": "#004488"}  # Paul Tol high-contrast
MARKER = {"interFoam": "o", "interIsoFoam": "s", "interIsoGoam": "D"}
CO = [0.5, 0.2, 0.1, 0.05]
TIMES = r"$\times$"

data = {}
for line in open(os.path.join(HERE, "ddtCorr.txt")):
    w = line.split()
    if len(w) < 3 or w[2] not in ("dHfit", "r0"):
        continue
    v = [float(x) for x in w[3:] if x != "/"]
    data[(w[0], w[1], w[2])] = (v[0::2], v[1::2])

fig, axes = plt.subplots(1, 3, figsize=(7.0, 3.0))
for ax, mesh, lab in zip(axes[:2], ["190x40", "380x80"], ["(a)", "(b)"]):
    for app in COLOUR:
        on, off = data[(mesh, app, "dHfit")]
        ax.plot(CO, on, "-", color=COLOUR[app], marker=MARKER[app])
        ax.plot(CO, off, "--", color=COLOUR[app], marker=MARKER[app], mfc="white")
    ax.axhline(0, color="#1f1f1e", lw=0.6)
    ax.set_title(f"{lab} {mesh.replace('x', TIMES)}", loc="left")
    ax.set_ylabel(r"$\Delta H$ after 15 periods [%]")
    ax.set_ylim(-36, 14)
ax = axes[2]
for mesh, c, m in [("190x40", "#6699CC", "o"), ("380x80", COLOUR["interIsoGoam"], "D")]:
    on, off = data[(mesh, "interIsoGoam", "r0")]
    ax.plot(CO, on, "-", color=c, marker=m, label=mesh.replace('x', TIMES))
    ax.plot(CO, off, "--", color=c, marker=m, mfc="white")
ax.axhline(1, color="#1f1f1e", lw=0.6)
ax.set_title("(c) interIsoGoam", loc="left")
ax.set_ylabel(r"$u/u_{\mathrm{exact}}$, top water cell")
ax.legend(loc="upper left", handlelength=1.6)
for ax in axes:
    ax.set_xscale("log")
    ax.set_xticks(CO)
    ax.set_xticklabels([f"{c:g}" for c in CO])
    ax.minorticks_off()
    ax.grid(True, lw=0.4, alpha=0.5)
# Solver and line-style legend below (a) and (b)
h = [plt.Line2D([], [], color=COLOUR[a], marker=MARKER[a], label=a) for a in COLOUR]
h += [plt.Line2D([], [], color="#6b6b6b", ls="-", label="ddtCorr on"),
      plt.Line2D([], [], color="#6b6b6b", ls="--", marker="o", mfc="white", label="ddtCorr off")]
fig.tight_layout(rect=(0, 0.15, 1, 1))
# Shared x label just below the tick labels, and the legend just below the label
r = fig.canvas.get_renderer(); toFig = fig.transFigure.inverted()
y0 = min(ax.get_tightbbox(r).transformed(toFig).y0 for ax in axes)
t = fig.text(0.5, y0 - 0.01, "maximum Courant number", ha="center", va="top")
fig.legend(handles=h, loc="upper center", ncol=5, handlelength=1.6, columnspacing=1.0,
           bbox_to_anchor=(0.5, t.get_window_extent(r).transformed(toFig).y0 - 0.01))
fig.savefig(os.path.join(HERE, "ddtCorr.pdf"))
fig.savefig(os.path.join(HERE, "ddtCorr.png"), dpi=150)
