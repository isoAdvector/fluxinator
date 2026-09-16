#!/usr/bin/env python
#
# Plots the maximum velocity magnitude against time for the tilted box
# (hydrostatic equilibrium, exact solution U = 0) for the gravity models
# compared in the paper. gravityRecon (an older version of geometric) is left
# out.
#
# Usage:
#   ./plotMaxU.py [--scanDir gravityModelScan_1000deltaT] [--out maxU_tiltedBox]

import argparse
import os

import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 10,
    "axes.linewidth": 0.8,
    "lines.linewidth": 1.4,
})

# (case, label, colour, line style) in the order of the paper figure
CASES = [
    ("0000", "interFoam, basic", "#2a78d6", "-"),
    ("0001", "interIsoFoam, basic", "#eb6834", "-"),
    ("0002", "interIsoGoam, constantHeight", "#555555", "--"),
    ("0003", "interIsoGoam, geometric", "#1baf7a", "-"),
]

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--scanDir", default="gravityModelScan_1000deltaT")
parser.add_argument("--out", default="maxU_tiltedBox", help="Output file name without extension")
args = parser.parse_args()

fig, ax = plt.subplots(figsize=(5.0, 3.2))
for case, label, colour, style in CASES:
    f = os.path.join(args.scanDir, case, "postProcessing", "fieldMinMax1", "0", "fieldMinMax.dat")
    t, umax = np.loadtxt(f, comments="#", unpack=True, usecols=(0, 2), ndmin=2)
    ax.semilogy(t, umax, style, color=colour, label=label)
    print(f"{case} {label:30s} first {umax[0]:.3e}  end {umax[-1]:.3e}  max {umax.max():.3e}")

ax.set_xlabel(r"$t$ [s]")
ax.set_ylabel(r"max$|\mathbf{U}|$ [m/s]")
ax.set_xlim(0, t[-1])
ax.grid(True, which="major", linewidth=0.4, alpha=0.6)
ax.legend(loc="center right", fontsize=8, frameon=False)
fig.tight_layout()
fig.savefig(args.out + ".pdf")
fig.savefig(args.out + ".png", dpi=200)
print("Wrote", args.out + ".pdf and .png")
