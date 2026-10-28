"""Schematic of the height in the gravity term on a tilted mesh with a flat interface.
(a) basic model: the height of the face centres of the faces across which the density changes;
(b) interface height: the centroids of the interface segments, all at the interface height.
Cells are square, rotated by 20 degrees (exaggerated), cut by the horizontal plane z = 0.
Writes figures/height.pdf. Usage: python3 figHeight.py"""
import numpy as np
from matplotlib.patches import Polygon
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({"font.family": "serif", "font.size": 9})
th = np.radians(20); R = np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
N, h, zw = 6, 1.0, 0.13
cells = {}
for i in range(N):
    for j in range(N):
        c = R @ (np.array([[i, j], [i + 1, j], [i + 1, j + 1], [i, j + 1]]).T - N/2)*h
        cells[i, j] = c.T


def clip_below(P, z0):
    out = []
    for k in range(len(P)):
        a, b = P[k], P[(k + 1) % len(P)]
        if a[1] <= z0: out.append(a)
        if (a[1] - z0)*(b[1] - z0) < 0:
            s = (z0 - a[1])/(b[1] - a[1]); out.append(a + s*(b - a))
    return np.array(out)


def area(P):
    x, y = P[:, 0], P[:, 1]
    return 0.5*abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))


alpha = {k: (area(clip_below(P, zw))/area(P) if len(clip_below(P, zw)) > 2 else 0.0) for k, P in cells.items()}
fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.4))
for ax, title in zip(axes, ["(a) basic model: face-centre heights",
                            "(b) geometric model: interface heights"]):
    for k, P in cells.items():
        ax.add_patch(Polygon(P, fc="#f3f1ec", ec="#9a9893", lw=0.6))
        W = clip_below(P, zw)
        if len(W) > 2: ax.add_patch(Polygon(W, fc="#9cc3ea", ec="none"))
    ax.axhline(zw, color="#2a78d6", lw=1.2)
    for (i, j), P in cells.items():
        for di, dj in [(1, 0), (0, 1)]:
            nb = (i + di, j + dj)
            if nb not in cells or abs(alpha[i, j] - alpha[nb]) < 1e-9: continue
            # shared face of (i,j) and nb
            f = (P[1], P[2]) if di == 1 else (P[3], P[2])
            fc = 0.5*(f[0] + f[1])
            if "basic" in title:
                ax.plot(*fc, "o", color="#eb6834", ms=5)
                ax.plot([fc[0] - 0.25, fc[0] + 0.25], [fc[1], fc[1]], color="#eb6834", lw=0.8)
    if "geometric" in title:
        for k, P in cells.items():
            if 0 < alpha[k] < 1:
                pts = []
                for m in range(4):
                    a, b = P[m], P[(m + 1) % 4]
                    if (a[1] - zw)*(b[1] - zw) < 0:
                        s = (zw - a[1])/(b[1] - a[1]); pts.append(a + s*(b - a))
                if len(pts) < 2: continue
                cpt = np.mean(pts, axis=0)
                ax.plot(*cpt, "D", color="#1baf7a", ms=5)
    ax.set_xlim(-3.0, 3.0); ax.set_ylim(-1.9, 2.1); ax.set_aspect("equal"); ax.axis("off")
    ax.set_title(title, fontsize=9, loc="left")
axes[0].text(-2.9, -1.8, "water", color="#2a78d6"); axes[0].text(-2.9, 1.85, "air", color="#6b6b6b")
save = lambda n: (fig.savefig(f"figures/{n}.pdf", bbox_inches="tight"), fig.savefig(f"figures/{n}.png", dpi=150, bbox_inches="tight"))
save("height")
