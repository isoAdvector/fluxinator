#!/usr/bin/env python
#
# Plots the maximum velocity magnitude after each of the 10 time steps for the
# hydrostatic interface on non-orthogonal meshes: (a) the skewed box and (b) the
# tilted box with a refined region, all with corrected laplacian and snGrad
# schemes unless stated. Colour identifies the gravity model, line style the
# start: zero initial p_rgh without non-orthogonal correctors, zero initial p_rgh
# with correctors, hydrostatic initial p_rgh, and uncorrected schemes.
#
# Usage (from the test directory, after the scans in skewBox and
# tiltedBoxRefined have been run):
#   ./plotNonOrthMaxU.py [--out nonOrthMaxU]

import argparse
import os

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 10,
    "axes.linewidth": 0.8,
    "lines.linewidth": 1.4,
})

HERE = os.path.dirname(os.path.abspath(__file__))

MODELS = {
    "basic": ("interIsoFoam, basic", "#eb6834"),
    "constantHeight": ("interIsoGoam, constantHeight", "#555555"),
    "geometric1": ("interIsoGoam, geometric, nLayers 1", "#1baf7a"),
    "geometric2": ("interIsoGoam, geometric, nLayers 2", "#7b4fd6"),
}

STARTS = {
    "zero": ("zero $p_{rgh}$, no correctors", "-", "o"),
    "correctors": ("zero $p_{rgh}$, 2 correctors", "--", "s"),
    "hydrostatic": ("hydrostatic $p_{rgh}$, no correctors", ":", "^"),
    "uncorrected": ("uncorrected schemes", "-.", "D"),
}

# (case directory, model, start)
SKEWBOX = [
    ("skewBox/nonOrthCorrScan/0005", "basic", "zero"),
    ("skewBox/nonOrthCorrScan/0007", "constantHeight", "zero"),
    ("skewBox/nonOrthCorrScan/0008", "constantHeight", "correctors"),
    ("skewBox/hydrostaticInitScan/0002", "constantHeight", "hydrostatic"),
    ("skewBox/nonOrthCorrScan/0000", "geometric1", "zero"),
    ("skewBox/nonOrthCorrScan/0002", "geometric1", "correctors"),
    ("skewBox/hydrostaticInitScan/0000", "geometric1", "hydrostatic"),
    ("skewBox/skewBoxScan_z/0005", "geometric1", "uncorrected"),
]

REFINED = [
    ("tiltedBoxRefined/zeroStartScan/0006", "basic", "zero"),
    ("tiltedBoxRefined/zeroStartScan/0004", "constantHeight", "zero"),
    ("tiltedBoxRefined/zeroStartScan/0008", "constantHeight", "correctors"),
    ("tiltedBoxRefined/hydrostaticInitScan/0002", "constantHeight", "hydrostatic"),
    ("tiltedBoxRefined/zeroStartScan/0000", "geometric1", "zero"),
    ("tiltedBoxRefined/zeroStartScan/0002", "geometric1", "correctors"),
    ("tiltedBoxRefined/hydrostaticInitScan/0000", "geometric1", "hydrostatic"),
    ("tiltedBoxRefined/zeroStartScan/0007", "geometric1", "uncorrected"),
    ("tiltedBoxRefined/nLayersZeroStartScan/0002", "geometric2", "zero"),
    ("tiltedBoxRefined/nLayersZeroStartScan/0001", "geometric2", "correctors"),
    ("tiltedBoxRefined/nLayersHydrostaticScan/0001", "geometric2", "hydrostatic"),
]

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--out", default="nonOrthMaxU", help="Output file name without extension")
args = parser.parse_args()

fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.4), sharey=True)
for ax, cases, title in zip(axes, [SKEWBOX, REFINED],
                            ["(a) skewed box", "(b) tilted box with refined region"]):
    print(title)
    for case, model, start in cases:
        f = os.path.join(HERE, case, "postProcessing", "fieldMinMax1", "0", "fieldMinMax.dat")
        t, umax = np.loadtxt(f, comments="#", unpack=True, usecols=(0, 2), ndmin=2)
        steps = np.rint(t/0.01).astype(int)
        _, colour = MODELS[model]
        _, style, marker = STARTS[start]
        ax.semilogy(steps, umax, style, color=colour, marker=marker, markersize=3.5)
        print(f"  {case:45s} {model:15s} {start:12s} step 1 {umax[0]:.2e}  step 10 {umax[-1]:.2e}")
    ax.set_title(title, fontsize=10)
    ax.set_xlabel("time step")
    ax.set_xlim(0.5, 10.5)
    ax.set_xticks(range(1, 11))
    ax.grid(True, which="major", linewidth=0.4, alpha=0.6)
axes[0].set_ylabel(r"max$|\mathbf{U}|$ [m/s]")
axes[0].set_ylim(1e-12, 1e2)

modelHandles = [Line2D([], [], color=c, label=l) for l, c in MODELS.values()]
startHandles = [Line2D([], [], color="black", linestyle=s, marker=m, markersize=3.5, label=l)
                for l, s, m in STARTS.values()]
fig.legend(handles=modelHandles, loc="lower left", bbox_to_anchor=(0.08, 0.0), fontsize=8,
           frameon=False, ncol=1)
fig.legend(handles=startHandles, loc="lower left", bbox_to_anchor=(0.53, 0.0), fontsize=8,
           frameon=False, ncol=1)
fig.tight_layout(rect=(0, 0.25, 1, 1))
fig.savefig(os.path.join(HERE, args.out + ".pdf"))
fig.savefig(os.path.join(HERE, args.out + ".png"), dpi=200)
print("Wrote", args.out + ".pdf and .png")
