#!/usr/bin/env python
#
# Plots the maximum velocity magnitude against time for the sheared box
# (hydrostatic equilibrium, exact solution U = 0), one panel per scan, in the
# same style as test/tiltedBox/plotMaxU.py. gravityRecon (an older version of
# geometric) is left out. The geometric model with uncorrected (non-orthogonal
# correction switched off) laplacian and snGrad schemes is added as 0003_uncor.
#
# Usage:
#   ./plotMaxU.py [--out maxU_shearedBox]

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

SCANS = [
    ("gravityModelScan_10deltaT_x", "(a) mesh rotated 90°: rows level, columns sheared"),
    ("gravityModelScan_10deltaT_z", "(b) mesh not rotated: rows sheared, columns vertical"),
]

# (case, label, colour, line style, marker)
CASES = [
    ("0000", "interFoam, basic", "#2a78d6", "-", "o"),
    ("0001", "interIsoFoam, basic", "#eb6834", "-", "s"),
    ("0002", "interIsoGoam, constantHeight", "#555555", "--", "^"),
    ("0003", "interIsoGoam, geometric", "#1baf7a", "-", "D"),
    ("0003_uncor", "interIsoGoam, geometric, uncorrected", "#1baf7a", ":", "d"),
]

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--out", default="maxU_shearedBox", help="Output file name without extension")
args = parser.parse_args()

fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.4), sharey=True)
for ax, (scan, title) in zip(axes, SCANS):
    print(scan)
    for case, label, colour, style, marker in CASES:
        f = os.path.join(scan, case, "postProcessing", "fieldMinMax1", "0", "fieldMinMax.dat")
        if not os.path.isfile(f):
            print(f"  {case}: missing {f}")
            continue
        t, umax = np.loadtxt(f, comments="#", unpack=True, usecols=(0, 2), ndmin=2)
        ax.semilogy(t, umax, style, color=colour, marker=marker, markersize=4,
                    mfc="white" if style == ":" else colour, label=label)
        print(f"  {case:10s} {label:38s} first {umax[0]:.3e}  end {umax[-1]:.3e}  max {umax.max():.3e}")
    ax.set_title(title, fontsize=9, loc="left")
    ax.set_xlabel(r"$t$ [s]")
    ax.set_xlim(0, 0.1)
    ax.grid(True, which="major", linewidth=0.4, alpha=0.6)
axes[0].set_ylabel(r"max$|\mathbf{U}|$ [m/s]")
axes[0].set_ylim(1e-11, 1e2)
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc="upper center", ncol=3, fontsize=8, frameon=False,
           bbox_to_anchor=(0.5, 1.0))
fig.tight_layout(rect=(0, 0, 1, 0.86))
fig.savefig(args.out + ".pdf")
fig.savefig(args.out + ".png", dpi=200)
print("Wrote", args.out + ".pdf and .png")
