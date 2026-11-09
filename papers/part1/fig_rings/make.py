# Copyright (C) 2026 Johan Roenby
# SPDX-License-Identifier: GPL-3.0-or-later
"""Sketch of the rings of the geometric gravity model on a square mesh with a straight,
tilted interface: interface cells, ring 1 (point neighbours of interface cells) and ring 2
(point neighbours of ring 1), coloured as in Figure fig:reach; diamonds mark the centroids of the
interface segments. Pure geometry, no run data.
Part I, Figure fig:rings. Output: rings.pdf (and .png). Usage: python3 make.py"""
import os
import numpy as np
from matplotlib.patches import Patch, Rectangle
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))

plt.rcParams.update({"font.family": "serif", "font.size": 8.5})
# Ring colours of fig_reach/make.py
COLOURS = {0: (0.95, 0.76, 0.31), 1: (0.60, 0.82, 0.71), 2: (0.84, 0.15, 0.16)}
LABELS = {0: "interface cells", 1: "ring 1", 2: "ring 2"}
LINE = "#8a8883"

NX, NZ = 10, 9
Z0, SLOPE = 3.3, 0.22  # interface z = Z0 + SLOPE*x, in cell widths


def eta(x):
    return Z0 + SLOPE * x


# Interface cells: the line crosses the interior of the cell
ring = -np.ones((NX, NZ), dtype=int)
for i in range(NX):
    lo, hi = sorted((eta(i), eta(i + 1)))
    for k in range(NZ):
        if lo < k + 1 and hi > k:
            ring[i, k] = 0

# Ring n: point neighbours of ring n-1 not in a lower ring
for n in (1, 2):
    prev = np.argwhere(ring == n - 1)
    for i, k in prev:
        for di in (-1, 0, 1):
            for dk in (-1, 0, 1):
                a, b = i + di, k + dk
                if 0 <= a < NX and 0 <= b < NZ and ring[a, b] < 0:
                    ring[a, b] = n

fig, ax = plt.subplots(figsize=(3.0, 2.4))
for i in range(NX):
    for k in range(NZ):
        r = ring[i, k]
        fc = COLOURS.get(r, "white")
        ax.add_patch(Rectangle((i, k), 1, 1, fc=fc, ec=LINE, lw=0.4))
for i, k in np.argwhere(ring == 2):
    ax.add_patch(Rectangle((i, k), 1, 1, fc="none", ec="black", lw=0.9, zorder=3))

# Segment centroids: midpoint of the part of the line inside each interface cell
for i, k in np.argwhere(ring == 0):
    xa = max(i, (k - Z0) / SLOPE)
    xb = min(i + 1, (k + 1 - Z0) / SLOPE)
    xm = 0.5 * (xa + xb)
    ax.plot(xm, eta(xm), "D", color="black", ms=2.2, zorder=5)

x = np.array([0, NX])
ax.plot(x, eta(x), color="#2a78d6", lw=1.2, zorder=4)

ax.set_xlim(0, NX)
ax.set_ylim(0, NZ)
ax.set_aspect("equal")
ax.axis("off")
handles = [Patch(fc=COLOURS[n], ec="black" if n == 2 else LINE, lw=0.6, label=LABELS[n])
           for n in (0, 1, 2)]
ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.0), ncol=3,
          frameon=False, handlelength=1.2, columnspacing=1.0, fontsize=7.5)

out = os.path.join(HERE, "rings")
fig.savefig(out + ".pdf", bbox_inches="tight", pad_inches=0.01)
fig.savefig(out + ".png", dpi=200, bbox_inches="tight", pad_inches=0.01)
print(f"Wrote {out}.pdf and .png")
