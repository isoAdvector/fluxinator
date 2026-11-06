#!/usr/bin/env python3
# Copyright (C) 2026 Johan Roenby
# SPDX-License-Identifier: GPL-3.0-or-later
#
# Maximum velocity magnitude after the first time step in the untilted box with a
# refined region, against the still water level z_w, which is scanned over one
# coarse cell (0.48 <= z_w <= 0.52).
#   (a) the cells at the right refinement boundary, the scanned band of levels and,
#       for z_w = 0.49, the cell where the standard gravity term gives its largest
#       velocity;
#   (b) uncorrected Laplacian and face-normal gradient schemes, zero initial p_rgh;
#   (c) corrected schemes, hydrostatic initial p_rgh.
# interFoam and interIsoFoam are identical in the first step (the interface has not
# moved), so one line shows the basic gravity term of both.
#
# Usage:
# Part I, Figure fig:levelScan. Input: run/refinedBox/waterLevel_10steps.
# Output: levelScan.pdf. Usage:
#   python3 make.py

import argparse
import glob
import os
import re
import sys

import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
from common import case  # noqa: E402
import numpy as np
from matplotlib.patches import Rectangle

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 8,
    "axes.linewidth": 0.7,
    "lines.linewidth": 1.2,
})

# Colours of the paper's tilted box figures
BASIC = ("interFoam / interIsoFoam, basic", "#2a78d6", "-", "o")
CONST = ("interIsoGoam, constantHeight", "#555555", "--", None)
GEO1 = ("interIsoGoam, geometric, $n_L = 1$", "#eda100", "-", "s")
GEO2 = ("interIsoGoam, geometric, $n_L = 2$", "#1baf7a", "-", "^")
PANELS = [
    ("(b) uncorrected schemes", [("iF_unc", BASIC), ("g1_unc", GEO1), ("cH_unc", CONST)]),
    ("(c) corrected, hydrostatic start",
     [("iF_hyd", BASIC), ("g1_hyd", GEO1), ("g2_hyd", GEO2), ("cH_hyd", CONST)]),
]

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--scanDir", default=case("refinedBox", "waterLevel_10steps"))
parser.add_argument("--out", default=os.path.join(HERE, "levelScan"), help="Output file name without extension")
args = parser.parse_args()


def firstStep(config):
    """(z_w, max|U| after the first step) over the scanned levels"""
    pts = []
    for case in glob.glob(os.path.join(args.scanDir, config + "_z*")):
        zw = float(re.search(r"_z([0-9.]+)$", case).group(1))
        fmm = os.path.join(case, "postProcessing", "fieldMinMax1", "0", "fieldMinMax.dat")
        _, umax = np.loadtxt(fmm, comments="#", unpack=True, usecols=(0, 2), ndmin=2)
        pts.append((zw, umax[0]))
    return np.array(sorted(pts))


fig, axes = plt.subplots(1, 3, figsize=(6.3, 2.7), gridspec_kw={"width_ratios": [0.9, 1, 1]})

# (a) cells at the right refinement boundary x = 0.76
ax = axes[0]
zw = 0.49
ax.add_patch(Rectangle((0.70, 0.40), 0.14, zw - 0.40, color="#cfe3f7", lw=0))
for z0 in (0.44, 0.48, 0.52, 0.56):
    ax.plot([0.76, 0.84], [z0, z0], color="k", lw=0.7)
for x0 in (0.76, 0.80, 0.84):
    ax.plot([x0, x0], [0.40, 0.60], color="k", lw=0.7)
for z0 in np.arange(0.40, 0.601, 0.02):
    ax.plot([0.70, 0.76], [z0, z0], color="k", lw=0.5)
for x0 in (0.70, 0.72, 0.74):
    ax.plot([x0, x0], [0.40, 0.60], color="k", lw=0.5)
ax.plot([0.70, 0.84], [zw, zw], color="#2a78d6", lw=1.4)
ax.add_patch(Rectangle((0.74, 0.50), 0.02, 0.02, fill=False, edgecolor="#d62728", lw=1.8, zorder=5))
ax.annotate("", xy=(0.855, 0.48), xytext=(0.855, 0.52),
            arrowprops=dict(arrowstyle="<->", lw=0.8, shrinkA=0, shrinkB=0))
ax.text(0.862, 0.50, "$z_w$", va="center", fontsize=8)
ax.text(0.73, 0.585, "fine", ha="center", fontsize=7, bbox=dict(fc="white", ec="none", pad=0.5))
ax.text(0.80, 0.585, "coarse", ha="center", fontsize=7, bbox=dict(fc="white", ec="none", pad=0.5))
ax.set_xlim(0.695, 0.885)
ax.set_ylim(0.40, 0.60)
ax.set_aspect("equal")
ax.set_xticks([0.72, 0.76, 0.80])
ax.tick_params(labelsize=7)
ax.set_yticks([0.44, 0.48, 0.52, 0.56])
ax.set_xlabel("$x$ [m]")
ax.set_ylabel("$z$ [m]")
ax.set_title("(a) refinement boundary", fontsize=8)

for ax, (title, lines) in zip(axes[1:], PANELS):
    for config, (label, colour, style, marker) in lines:
        pts = firstStep(config)
        if pts.size == 0:
            continue
        ax.semilogy(pts[:, 0], pts[:, 1], style, color=colour, marker=marker, ms=3, label=label)
    ax.set_xlim(0.479, 0.521)
    ax.set_ylim(1e-14, 30)
    ax.set_xticks([0.48, 0.49, 0.50, 0.51, 0.52])
    ax.set_xlabel("water level $z_w$ [m]")
    ax.set_title(title, fontsize=8)
    ax.grid(True, which="major", lw=0.3, color="#cccccc")
axes[1].set_ylabel(r"$\max|\mathbf{U}|$ after one step [m/s]")
axes[2].set_yticklabels([])

handles, labels = [], []
for ax in axes[1:]:
    for h, l in zip(*ax.get_legend_handles_labels()):
        if l not in labels:
            handles.append(h)
            labels.append(l)
fig.legend(handles, labels, loc="lower center", ncol=2, frameon=False, fontsize=7,
           bbox_to_anchor=(0.6, -0.02))
fig.tight_layout(rect=(0, 0.12, 1, 1))
fig.savefig(args.out + ".pdf")
print("Wrote", args.out + ".pdf")
