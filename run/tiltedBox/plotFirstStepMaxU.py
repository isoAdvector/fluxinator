#!/usr/bin/env python
#
# Plots the maximum velocity magnitude after the first time step against the
# time step for interFoam and interIsoFoam (interIsoGoam with the basic gravity
# model) in the tilted box, from the one-step scan deltaTScan (deltaTScan.csv).
# A line of slope one is added for reference.
#
# Usage:
#   ./plotFirstStepMaxU.py [--scanDir deltaTScan] [--out firstStepMaxU_tiltedBox]

import argparse
import os
import re

import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 10,
    "axes.linewidth": 0.8,
    "lines.linewidth": 1.4,
})

STYLE = {
    # interFoam drawn first with larger filled markers so that it stays visible
    # beneath the identical interIsoFoam values
    "interFoam": ("interFoam, basic", "#2a78d6", "o", "-", 10),
    "interIsoGoam": ("interIsoFoam, basic", "#eb6834", "s", "--", 5),
}

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--scanDir", default="deltaTScan")
parser.add_argument("--out", default="firstStepMaxU_tiltedBox", help="Output file name without extension")
args = parser.parse_args()

# Case parameters from the scan's parmScannerParms
data = {}
with open(os.path.join(args.scanDir, "parmScannerParms")) as f:
    for line in f:
        case = line.split(":")[0].split("/")[-1]
        parms = dict(re.findall(r"(\w+):\s*(\S+)", line.split(":", 1)[1]))
        fmm = os.path.join(args.scanDir, case, "postProcessing", "fieldMinMax1", "0", "fieldMinMax.dat")
        t, umax = np.loadtxt(fmm, comments="#", unpack=True, usecols=(0, 2), ndmin=2)
        data.setdefault(parms["application"], []).append((float(parms["deltaT"]), float(umax[0])))

fig, ax = plt.subplots(figsize=(4.2, 3.2))
for app, (label, colour, marker, style, size) in STYLE.items():
    pts = sorted(data.get(app, []))
    dt = np.array([p[0] for p in pts])
    u = np.array([p[1] for p in pts])
    ax.loglog(dt, u, style, color=colour, marker=marker, mfc="white" if style == "--" else colour,
              markersize=size, label=label)
    print(app)
    for d, v in pts:
        print(f"  dt={d:8.1e}  max|U|={v:.6e}  max|U|/dt={v/d:.6e}")
# Slope-one reference through the interFoam point at dt = 1e-2
ref = sorted(data["interFoam"])
d0, u0 = min(ref, key=lambda p: abs(np.log10(p[0]) + 2))
dd = np.array([1e-4, 1e0])
ax.loglog(dd, u0*dd/d0, ":", color="k", lw=1.0, label=r"$\propto \Delta t$")
ax.set_xlabel(r"$\Delta t$ [s]")
ax.set_ylabel(r"max$|\mathbf{U}|$ after first time step [m/s]")
ax.grid(True, which="major", linewidth=0.4, alpha=0.6)
ax.legend(loc="upper left", fontsize=8, frameon=False)
fig.tight_layout()
fig.savefig(args.out + ".pdf")
fig.savefig(args.out + ".png", dpi=200)
print("Wrote", args.out + ".pdf and .png")
