# Copyright (C) 2026 Johan Roenby
# SPDX-License-Identifier: GPL-3.0-or-later
"""Sketch of the layers of the geometric gravity model on a square mesh with a straight,
tilted interface: interface cells, layer 1 (point neighbours of interface cells) and layer 2
(point neighbours of layer 1), coloured as in Figure fig:reach; dots mark the centroids of the
interface segments. Pure geometry, no run data.
Part I, Figure fig:layers. Output: layers.pdf (and .png). Usage: python3 make.py"""
import os
import sys
import numpy as np
from matplotlib.patches import Patch, Rectangle
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))

sys.path.insert(0, os.path.join(HERE, ".."))
from common import paper_style  # noqa: E402

paper_style()
# Layer colours of fig_reach/make.py
COLOURS = {0: (0.95, 0.76, 0.31), 1: (0.60, 0.82, 0.71), 2: (0.84, 0.15, 0.16)}
LABELS = {0: "interface cells", 1: "layer 1", 2: "layer 2"}
LINE = "#8a8883"

NX, NZ = 10, 9
Z0, SLOPE = 3.3, 0.22  # interface z = Z0 + SLOPE*x, in cell widths


def eta(x):
    return Z0 + SLOPE * x


# Interface cells: the line crosses the interior of the cell
layer = -np.ones((NX, NZ), dtype=int)
for i in range(NX):
    lo, hi = sorted((eta(i), eta(i + 1)))
    for k in range(NZ):
        if lo < k + 1 and hi > k:
            layer[i, k] = 0

# Layer n: point neighbours of layer n-1 not in a lower layer
for n in (1, 2):
    prev = np.argwhere(layer == n - 1)
    for i, k in prev:
        for di in (-1, 0, 1):
            for dk in (-1, 0, 1):
                a, b = i + di, k + dk
                if 0 <= a < NX and 0 <= b < NZ and layer[a, b] < 0:
                    layer[a, b] = n

fig, ax = plt.subplots(figsize=(2.9, 2.75))  # 0.45 of the text width after the tight crop
for i in range(NX):
    for k in range(NZ):
        r = layer[i, k]
        fc = COLOURS.get(r, "white")
        ax.add_patch(Rectangle((i, k), 1, 1, fc=fc, ec=LINE, lw=0.4))
for i, k in np.argwhere(layer == 2):
    ax.add_patch(Rectangle((i, k), 1, 1, fc="none", ec="black", lw=0.9, zorder=3))

# Segment centroids: midpoint of the part of the line inside each interface cell
for i, k in np.argwhere(layer == 0):
    xa = max(i, (k - Z0) / SLOPE)
    xb = min(i + 1, (k + 1 - Z0) / SLOPE)
    xm = 0.5 * (xa + xb)
    ax.plot(xm, eta(xm), "o", color="black", ms=2.0, zorder=5)

x = np.array([0, NX])
ax.plot(x, eta(x), color="#2a78d6", lw=1.2, zorder=4)

ax.set_xlim(0, NX)
ax.set_ylim(0, NZ)
ax.set_aspect("equal")
ax.axis("off")
handles = [Patch(fc=COLOURS[n], ec="black" if n == 2 else LINE, lw=0.6, label=LABELS[n])
           for n in (0, 1, 2)]
ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.0), ncol=3,
          frameon=False, handlelength=1.0, handleheight=1.0, handletextpad=0.4, columnspacing=0.8)

out = os.path.join(HERE, "layers")
fig.savefig(out + ".pdf", bbox_inches="tight", pad_inches=0.01)
fig.savefig(out + ".png", dpi=200, bbox_inches="tight", pad_inches=0.01)
print(f"Wrote {out}.pdf and .png")
