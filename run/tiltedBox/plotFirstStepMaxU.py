#!/usr/bin/env python
#
# Maximum velocity magnitude after the first time step in the tilted box, against
# the time step, for meshes from 25 x 25 to 200 x 200 cells and for both
# interFoam and interIsoFoam (run as interIsoGoam with the basic gravity model).
#
# Each mesh gives a straight line of slope one, so max|U|/dt is a spurious
# acceleration independent of the time step, and the lines of the four meshes
# lie almost on top of each other: the error of the standard gravity term does
# not fall with the cell size, as Section 2.3 predicts. The two solvers give
# identical values at the first step, before the interface has moved.
#
# Usage:
#   ./plotFirstStepMaxU.py [--scanDir meshDeltaTScan] [--out firstStepMaxU_tiltedBox]

import argparse
import os
import re

import matplotlib.pyplot as plt
import numpy as np

# Sized to sit beside maxU_tiltedBox.pdf at half the paper's text width
plt.rcParams.update({
    "font.family": "serif",
    "font.size": 8,
    "axes.linewidth": 0.7,
    "lines.linewidth": 1.2,
})

G = 9.81
# interFoam drawn first with larger filled markers so that it stays visible
# beneath the identical interIsoFoam values
SOLVER = {"interFoam":    ("interFoam", "o", "-", 7, True),
          "interIsoGoam": ("interIsoFoam", "s", "--", 3, False)}
MESH_COLOUR = {25: "#9dc3f0", 50: "#5f9fe0", 100: "#2a78d6", 200: "#12457f"}

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--scanDir", default="meshDeltaTScan")
parser.add_argument("--out", default="firstStepMaxU_tiltedBox", help="Output file name without extension")
args = parser.parse_args()

data = {}
with open(os.path.join(args.scanDir, "parmScannerParms")) as f:
    for line in f:
        case = line.split(":")[0].split("/")[-1]
        p = dict(re.findall(r"(\w+):\s*(\S+)", line.split(":", 1)[1]))
        fmm = os.path.join(args.scanDir, case, "postProcessing", "fieldMinMax1", "0",
                           "fieldMinMax.dat")
        _, umax = np.loadtxt(fmm, comments="#", unpack=True, usecols=(0, 2), ndmin=2)
        data.setdefault((int(p["nx"]), p["application"]), []).append(
            (float(p["deltaT"]), float(umax[0])))

fig, ax = plt.subplots(figsize=(3.05, 2.4))
for N in sorted(MESH_COLOUR):
    for app, (label, marker, style, size, fill) in SOLVER.items():
        pts = sorted(data.get((N, app), []))
        if not pts:
            continue
        dt = np.array([p[0] for p in pts])
        u = np.array([p[1] for p in pts])
        ax.loglog(dt, u, style, color=MESH_COLOUR[N], marker=marker,
                  mfc=MESH_COLOUR[N] if fill else "white", markersize=size,
                  mew=0.8, label=rf"${N}\times{N}$" if fill else None)

print("spurious acceleration max|U|/dt of the first time step")
for N in sorted(MESH_COLOUR):
    for app in SOLVER:
        a = np.array([u/dt for dt, u in data[(N, app)]])
        print(f"  {N:4d}^2  {app:13s} h = {1.0/N:7.4f} m  a = {a.mean():8.3f} m/s2 "
              f"= {a.mean()/G:.3f} g   (spread over the five dt: {a.ptp():.1e})")
    ai = np.array([u/dt for dt, u in data[(N, 'interFoam')]])
    ag = np.array([u/dt for dt, u in data[(N, 'interIsoGoam')]])
    print(f"        max relative difference between the solvers: "
          f"{np.abs(ai - ag).max()/ai.mean():.1e}")

mesh = [plt.Line2D([], [], color=MESH_COLOUR[N], lw=1.2, label=rf"${N}\times{N}$")
        for N in sorted(MESH_COLOUR)]
style = [plt.Line2D([], [], color="0.35", ls="-", marker="o", ms=5, label="interFoam"),
         plt.Line2D([], [], color="0.35", ls="--", marker="s", ms=3, mfc="white",
                    label="interIsoFoam")]
ax.legend(handles=mesh + style, loc="upper left", fontsize=6.5, frameon=False,
          ncol=2, columnspacing=0.8, handlelength=1.6, handletextpad=0.4,
          borderaxespad=0.3)
ax.set_xlabel(r"$\Delta t$ [s]")
ax.set_ylabel(r"max$|\mathbf{U}|$ after one time step [m/s]")
ax.grid(True, which="major", linewidth=0.4, alpha=0.6)
fig.tight_layout(pad=0.3)
fig.savefig(args.out + ".pdf")
fig.savefig(args.out + ".png", dpi=200)
print("Wrote", args.out + ".pdf and .png")
