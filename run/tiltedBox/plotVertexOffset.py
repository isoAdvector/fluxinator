#!/usr/bin/env python3
#
# Figure for the tilted box section: why particular cells are hit by the
# spurious currents of the standard gravity term, and which ones.
#
#   (a)-(c) the first-step velocity field of the interFoam runs at 5, 10 and 15
#           degrees, with the vertices the closed form points at.  The hot cells
#           move along the interface and get denser as the tilt grows, spaced by
#           cot(theta) columns;
#   (d)     the universal vertex loop imbalance Gamma(xi) - pure theory;
#   (e)     max|U| after the first time step against tilt angle, over the
#           interFoam runs of tiltScan_1deltaT, with the closed-form prediction.
#
# Everything drawn from a simulation comes from OpenFOAM runs.
#
# Usage:
#   ./plotVertexOffset.py [--out vertexOffset_tiltedBox]
#                        [--tiltScan tiltScan_1deltaT]
#                        [--fieldTilts 5 10 15]

import argparse

import matplotlib.pyplot as plt
import numpy as np

import tiltedBoxCase as tbc
from vertexOffsetTheory import (Gamma, GammaPeak, predictFromMesh, rankVertices,
                               C0, DT)

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 10,
    "axes.linewidth": 0.8,
    "lines.linewidth": 1.4,
})

BLUE, ORANGE = "#2a78d6", "#eb6834"
RAMP = ["#9dc3f0", "#3d86d8", "#12457f"]          # one hue, light -> dark
GAMMA_ANGLES = [5.0, 10.0, 15.0]
N, L = 25, 1.0
H = L/N

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--out", default="vertexOffset_tiltedBox")
ap.add_argument("--tiltScan", default="tiltScan_1deltaT")
ap.add_argument("--fieldTilts", type=float, nargs=3, default=[5.0, 10.0, 15.0])
args = ap.parse_args()

# Every case of the tilt scan, read once
runs = {}
for case in tbc.cases(args.tiltScan):
    r = tbc.firstStep(case, N)
    runs[round(r["theta"], 4)] = r

fig = plt.figure(figsize=(7.2, 5.2))
outer = fig.add_gridspec(2, 1, height_ratios=[1.0, 1.55], hspace=0.34)
top = outer[0].subgridspec(1, 3, wspace=0.12)
bot = outer[1].subgridspec(1, 2, wspace=0.34)
axF = [fig.add_subplot(top[0, j]) for j in range(3)]
axG = fig.add_subplot(bot[0, 0])
axU = fig.add_subplot(bot[0, 1])

# ------------------------------------------- (a)-(c) the velocity fields
for n, (ax, th) in enumerate(zip(axF, args.fieldTilts)):
    r = runs[round(th, 4)]
    V, mU = r["V"], r["magU"]
    ax.pcolormesh(V[..., 0], V[..., 1], mU/mU.max(), cmap="Blues", vmin=0, vmax=1,
                  shading="flat", rasterized=True, edgecolors="0.90", linewidth=0.10)
    ax.contour(r["C"][..., 0], r["C"][..., 1], r["alpha"], levels=[0.5],
               colors="0.25", linewidths=0.8)
    rk = rankVertices(V, r["zw"], r["theta"])
    k, i = np.where(rk["window"])
    ax.plot(V[k+1, i+1, 0], V[k+1, i+1, 1], "o", mfc="none", mec=ORANGE,
            mew=1.2, ms=7, ls="none")
    ax.set_aspect("equal")
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(r["zw"] - 6.5*H, r["zw"] + 6.5*H)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(rf"({'abc'[n]}) ${th:.0f}^\circ$,  "
                 rf"max$|\mathbf{{U}}|={mU.max():.2f}$ m/s", fontsize=8)
axF[0].plot([], [], "o", mfc="none", mec=ORANGE, mew=1.2, ms=7,
            label="vertex in the light-fluid window (16)")
axF[0].legend(fontsize=7.5, frameon=False, loc="upper center",
              bbox_to_anchor=(1.60, -0.02), handletextpad=0.3)

# ------------------------------------------------------- (d) Gamma(xi)
xi = np.linspace(-1.4, 1.4, 4001)*H
for th, col in zip(GAMMA_ANGLES, RAMP):
    axG.plot(xi/H, Gamma(xi, th, H), color=col, label=rf"${th:.0f}^\circ$")
t15 = np.radians(GAMMA_ANGLES[-1])
axG.axvspan(-np.cos(t15) - np.sin(t15), -np.cos(t15), color=ORANGE, alpha=0.16, lw=0)
axG.annotate("loop entirely\nin light fluid", xy=(-1.04, 6.0),
             xytext=(-1.34, 33.0), ha="left", va="top", fontsize=7.5,
             color="#8a3c14",
             arrowprops=dict(arrowstyle="->", color="#8a3c14", lw=0.7,
                             shrinkA=2, shrinkB=2))
for th, col in zip(GAMMA_ANGLES, RAMP):
    t = np.radians(th)
    axG.plot([0], [GammaPeak(th, H)], "o", color=col, ms=4.5, zorder=5)
    axG.plot([-np.cos(t)], [Gamma(-H*np.cos(t), th, H)[0]], "s", color=col,
             ms=4.5, mfc="white", zorder=5)
axG.axhline(0, color="0.4", lw=0.6)
axG.set_xlabel(r"vertex offset $\xi/h = (z_w - z_v)/h$")
axG.set_ylabel(r"loop imbalance $\Gamma$ [Pa]")
axG.set_xlim(-1.4, 1.4)
axG.set_ylim(-36, 42)
axG.grid(True, lw=0.4, alpha=0.5)
axG.legend(title="tilt", fontsize=8, title_fontsize=8, frameon=False,
           loc="upper right", handlelength=1.3, borderaxespad=0.2)
axG.set_title("(d) universal loop imbalance", fontsize=9)

# ---------------------------------------- (e) the tilt scan against theory
pts = np.array(sorted(
    (r["theta"], r["magU"].max(), predictFromMesh(r["V"], r["zw"], r["theta"])["maxU"])
    for r in runs.values()))
axU.plot(pts[:, 0], pts[:, 1], color=BLUE, lw=1.2, label=f"interFoam ({len(pts)} runs)")
axU.plot(pts[:, 0], pts[:, 2], color=ORANGE, lw=1.0, ls="--", label="closed form")
axU.plot(pts[:, 0], GammaPeak(pts[:, 0], H)*C0*DT/(H*1.0)/2, color="0.55", lw=0.9,
         ls=":", label=r"envelope $\propto \tan\theta$")
for th in args.fieldTilts:                       # tie (e) back to (a)-(c)
    axU.axvline(th, color="0.78", lw=0.7, ls="-", zorder=0)
for n, th in enumerate(args.fieldTilts):
    axU.annotate("abc"[n], xy=(th, 0.045), ha="center", va="center", fontsize=7.5,
                 color="0.35", bbox=dict(fc="white", ec="none", pad=0.6))
axU.set_xlabel(r"tilt angle $\theta$ [deg]")
axU.set_ylabel(r"max$|\mathbf{U}|$ [m/s]", labelpad=1)
axU.set_ylim(0, 1.45)
axU.grid(True, lw=0.4, alpha=0.5)
axU.legend(fontsize=7.5, frameon=False, loc="upper left", borderaxespad=0.3)
axU.set_title("(e) which tilt is unlucky", fontsize=9)

fig.savefig(args.out + ".pdf", bbox_inches="tight")
fig.savefig(args.out + ".png", dpi=200, bbox_inches="tight")
print("Wrote", args.out + ".pdf and .png")
