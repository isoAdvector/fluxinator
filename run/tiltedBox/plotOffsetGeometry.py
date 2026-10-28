#!/usr/bin/env python3
#
# Definition figure for the tilted box section: what the symbols of the offset
# argument mean geometrically.
#
#   (a) one cell of the tilted mesh - the two corner heights a and b above the
#       cell centre, the depth u = z_w - z_C, and the wetted fraction A(u) that
#       enters the universal loop imbalance;
#   (b) the four-cell loop at a vertex - the offset xi = z_w - z_v, and the band
#       of interface heights for which the loop lies entirely in the light fluid
#       while the imbalance is still non-zero;
#   (c) the fractional offset delta of each vertex column, a sawtooth of slope
#       tan(theta) that enters the light-fluid band every cot(theta) columns.
#
# Pure geometry: nothing here is read from a simulation.  The tilt is
# exaggerated in (a) and (b) so that a and b are distinguishable.
#
# Usage:
#   ./plotOffsetGeometry.py [--out offsetGeometry_tiltedBox] [--theta 25]

import argparse

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Polygon

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 10,
    "axes.linewidth": 0.8,
    "lines.linewidth": 1.4,
})

BLUE, ORANGE = "#2a78d6", "#eb6834"
WATER = "#cfe2f7"
INK = "0.25"

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--out", default="offsetGeometry_tiltedBox")
ap.add_argument("--theta", type=float, default=25.0, help="tilt drawn in (a),(b)")
ap.add_argument("--thetaC", type=float, default=5.0, help="tilt used in (c)")
ap.add_argument("-N", type=int, default=25)
args = ap.parse_args()

T = np.radians(args.theta)
C, S = np.cos(T), np.sin(T)
h = 1.0
E1 = np.array([C, -S])           # mesh axes, as in the case (rotate-y)
E3 = np.array([S,  C])
A_ = 0.5*h*(C - S)               # lower corner offset  a
B_ = 0.5*h*(C + S)               # upper corner offset  b


def cell(centre, hh=h):
    """Corners of a tilted square cell, anticlockwise."""
    o = np.array(centre, float)
    return np.array([o + s1*0.5*hh*E1 + s3*0.5*hh*E3
                     for s1, s3 in ((-1, -1), (1, -1), (1, 1), (-1, 1))])


def below(poly, y0):
    """The part of a convex polygon with y <= y0 (Sutherland-Hodgman)."""
    out = []
    n = len(poly)
    for i in range(n):
        p, q = poly[i], poly[(i+1) % n]
        pin, qin = p[1] <= y0, q[1] <= y0
        if pin:
            out.append(p)
        if pin != qin:
            t = (y0 - p[1])/(q[1] - p[1])
            out.append(p + t*(q - p))
    return np.array(out) if out else np.empty((0, 2))


fig = plt.figure(figsize=(7.2, 4.9))
outer = fig.add_gridspec(2, 1, height_ratios=[1.0, 0.86], hspace=0.16)
top = outer[0].subgridspec(1, 2, width_ratios=[0.86, 1.14], wspace=0.16)
axA = fig.add_subplot(top[0, 0])
axB = fig.add_subplot(top[0, 1])
axC = fig.add_subplot(outer[1])

# --------------------------------------------------- (a) one cell: a, b, A(u)
P = cell((0, 0))
zw = 0.09*h                                     # water plane, clear of +/-a
axA.add_patch(Polygon(below(P, zw), fc=WATER, ec="none", zorder=1))
axA.add_patch(Polygon(P, fc="none", ec=INK, lw=1.4, zorder=3))
axA.plot([-1.02, 1.10], [zw, zw], color=BLUE, lw=1.5, zorder=4)
axA.annotate(r"$z_w$", xy=(1.16, zw), va="center", fontsize=9, color=BLUE)
for y, lab in ((B_, "$+b$"), (A_, "$+a$"), (-A_, "$-a$"), (-B_, "$-b$")):
    axA.plot([-1.02, 1.10], [y, y], color="0.72", lw=0.6, ls=(0, (3, 2)), zorder=0)
    axA.annotate(lab, xy=(1.16, y), va="center", ha="left", fontsize=9, color="0.35")
axA.plot(P[:, 0], P[:, 1], "o", color="0.45", ms=3.2, ls="none", zorder=5)
axA.plot([0], [0], "o", color=INK, ms=4, zorder=5)
axA.annotate(r"$\mathbf{x}_C$", xy=(0.09, -0.20), fontsize=9)
axA.annotate("", xy=(-0.78, zw), xytext=(-0.78, 0),
             arrowprops=dict(arrowstyle="<->", color=ORANGE, lw=1.2))
axA.annotate(r"$u$", xy=(-0.86, 0.5*zw), ha="right", va="center",
             fontsize=9.5, color="#8a3c14")
axA.annotate(r"$\alpha = A(u)$", xy=(0.0, -1.06), ha="center",
             fontsize=9, color=BLUE)
axA.set_title("(a) one cell: the corner heights" "\n" r"$\pm a$, $\pm b$ and the cut $A(u)$",
              fontsize=9)
axA.set_xlim(-1.75, 1.75); axA.set_ylim(-1.25, 0.95)

# ------------------------------------- (b) the loop at a vertex: xi and window
V = np.array([0.0, 0.0])
for s1 in (-1, 1):
    for s3 in (-1, 1):
        axB.add_patch(Polygon(cell(V + s1*0.5*h*E1 + s3*0.5*h*E3),
                              fc="none", ec="0.80", lw=1.0, zorder=2))
low = cell(V + 0.5*h*E1 - 0.5*h*E3)             # the LOWEST cell: centre at -b
yBot, yTop = -(C + S)*h, -C*h                   # its two lowest corner heights
xL = -1.95
axB.axhspan(yBot, yTop, xmin=0.0, xmax=0.62, color=ORANGE, alpha=0.22, lw=0, zorder=0)
xi = 0.5*(yBot + yTop)
axB.add_patch(Polygon(below(low, xi), fc=WATER, ec="none", zorder=1))
axB.add_patch(Polygon(low, fc="none", ec=INK, lw=1.5, zorder=3))
axB.plot([xL, 0.95], [xi, xi], color=BLUE, lw=1.5, zorder=4)
axB.annotate(r"$z_w$", xy=(xL + 0.05, xi + 0.10), fontsize=9, color=BLUE)
axB.plot([0], [0], "o", color=ORANGE, ms=8, mfc="none", mew=1.8, zorder=6)
axB.annotate(r"vertex $v$", xy=(0.14, 0.10), fontsize=9, color="#8a3c14")
axB.annotate("", xy=(-0.62, xi), xytext=(-0.62, 0),
             arrowprops=dict(arrowstyle="<->", color=ORANGE, lw=1.3))
axB.annotate(r"$\xi$", xy=(-0.70, 0.5*xi), ha="right", va="center",
             fontsize=9.5, color="#8a3c14")
axB.plot([-0.62, 0], [0, 0], color=ORANGE, lw=0.6, ls=(0, (2, 2)), zorder=5)
axB.plot(low[:, 0], low[:, 1], "o", color="0.45", ms=3.2, ls="none", zorder=6)
for y in (yTop, yBot):
    axB.plot([xL, 1.26], [y, y], color="0.70", lw=0.6, ls=(0, (3, 2)), zorder=1)
axB.annotate(r"$-h\cos\theta$", xy=(2.00, yTop), va="center", fontsize=8, color="0.35")
axB.annotate(r"$-h(\cos\theta + \sin\theta)$", xy=(2.00, yBot), va="center",
             fontsize=8, color="0.35")
axB.annotate("", xy=(1.34, yBot), xytext=(1.34, yTop),
             arrowprops=dict(arrowstyle="<->", color="#8a3c14", lw=1.1))
axB.annotate(r"$h\sin\theta$", xy=(1.42, xi), ha="left", va="center", fontsize=8,
             color="#8a3c14")
axB.annotate(r"window (16): $\alpha \leq \frac{1}{2}\tan\theta$ in the lowest cell",
             xy=(0.10, yBot - 0.26), ha="center", va="top", fontsize=8, color=BLUE)
axB.set_title("(b) the four-cell loop at a vertex:" "\n"
              r"the offset $\xi$ and the window", fontsize=9)
axB.set_xlim(xL, 3.55); axB.set_ylim(-2.30, 1.30)

for ax in (axA, axB):
    ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)

# ------------------------------------------------ (c) the fractional-offset sawtooth
tC = np.radians(args.thetaC)
N = args.N
I = np.arange(1, N)
m = (I - N/2)*np.tan(tC) + N/2 + 1
delta = m - np.rint(m)
axC.axhspan(-np.tan(tC), 0.0, color=ORANGE, alpha=0.18, lw=0, zorder=0,
            label=r"light-fluid window (16)")
axC.axhline(0, color="0.5", lw=0.7, zorder=1)
axC.plot(I, delta, "o-", color=BLUE, ms=3.2, lw=0.9, zorder=3,
         label=r"$\delta(I)$, slope $\tan\theta$")
hot = (delta > -np.tan(tC)) & (delta <= 0)
axC.plot(I[hot], delta[hot], "o", color=ORANGE, ms=7, mfc="none", mew=1.5,
         zorder=4, label="hot columns")
if hot.sum() >= 2:
    x0, x1 = I[hot][0], I[hot][1]
    axC.annotate("", xy=(x0, 0.40), xytext=(x1, 0.40),
                 arrowprops=dict(arrowstyle="<->", color="0.35", lw=0.9))
    axC.annotate(r"$\cot\theta$ columns", xy=(0.5*(x0+x1), 0.43), ha="center",
                 va="bottom", fontsize=8, color="0.35")
axC.set_xlabel(r"vertex column $I$")
axC.set_ylabel(r"fractional offset $\delta$")
axC.set_ylim(-0.55, 0.58)
axC.set_xlim(0, N)
axC.grid(True, lw=0.4, alpha=0.5)
axC.legend(fontsize=7, frameon=False, loc="lower right", handletextpad=0.5,
           borderaxespad=0.3)
axC.set_title(rf"(c) offset lattice, $\theta={args.thetaC:.0f}^\circ$", fontsize=9)

fig.savefig(args.out + ".pdf", bbox_inches="tight")
fig.savefig(args.out + ".png", dpi=200, bbox_inches="tight")
print("Wrote", args.out + ".pdf and .png")
