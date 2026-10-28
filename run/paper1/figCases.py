"""Overview of the test cases of Part 1: mesh and initial interface (water shaded).
(a) tilted box, 25x25 rotated by 5 deg; (b) skewed box, 25x25 with rows rising 0.1 m over 1 m;
(c) box with a refined region (central 13x13 cells split 2x2), water level 0.49 m;
(d) standing wave, 50x50 with the same refinement, amplitude 0.05 m, wavelength 2 m;
(e) progressive wave in a periodic domain of one wavelength (mesh not drawn; 95x20 to 760x160).
Writes figures/cases.pdf. Usage: python3 figCases.py"""
import os
import sys
import numpy as np
from matplotlib.patches import Polygon, Rectangle
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "qwistWave"))
import kinematics as K  # noqa: E402

plt.rcParams.update({"font.family": "serif", "font.size": 8.5, "axes.titlesize": 8.5})
WATER, LINE = "#9cc3ea", "#8a8883"


def grid_lines(ax, X, Z, lw=0.25):
    for i in range(X.shape[0]): ax.plot(X[i], Z[i], color=LINE, lw=lw)
    for j in range(X.shape[1]): ax.plot(X[:, j], Z[:, j], color=LINE, lw=lw)


def water(ax, outline, eta, xr):
    dom = Polygon(outline, fc="none", ec="#1f1f1e", lw=0.8); ax.add_patch(dom)
    x = np.linspace(*xr, 400)
    poly = np.vstack([np.c_[x, eta(x)], [[xr[1], -10], [xr[0], -10]]])
    w = Polygon(poly, fc=WATER, ec="none", zorder=0); ax.add_patch(w); w.set_clip_path(dom)
    ln, = ax.plot(x, eta(x), color="#2a78d6", lw=1.0); ln.set_clip_path(dom)


fig = plt.figure(figsize=(7.2, 3.9))
gs = fig.add_gridspec(2, 4, height_ratios=[1.0, 0.55], hspace=0.25, wspace=0.12)
axes = [fig.add_subplot(gs[0, i]) for i in range(4)] + [fig.add_subplot(gs[1, :])]
# (a) tilted box
ax = axes[0]; th = np.radians(5); n = 25
s = np.linspace(-0.5, 0.5, n + 1); S, T = np.meshgrid(s, s)
X = 0.5 + S*np.cos(th) - T*np.sin(th); Z = 0.5 + S*np.sin(th) + T*np.cos(th)
grid_lines(ax, X, Z)
out = np.array([[X[0, 0], Z[0, 0]], [X[0, -1], Z[0, -1]], [X[-1, -1], Z[-1, -1]], [X[-1, 0], Z[-1, 0]]])
water(ax, out, lambda x: 0.5 + 0*x, (-0.2, 1.2)); ax.set_title("(a) tilted box")
# (b) skewed box
ax = axes[1]; x = np.linspace(0, 1, n + 1); z = np.linspace(0, 1, n + 1); X, Z = np.meshgrid(x, z)
Z = Z + 0.1*X - 0.05
grid_lines(ax, X, Z)
out = np.array([[0, -0.05], [1, 0.05], [1, 1.05], [0, 0.95]])
water(ax, out, lambda x: 0.5 + 0*x, (0, 1)); ax.set_title("(b) skewed box")


def refined(ax, n, eta, title):
    h = 1.0/n; lo, hi = 0.24, 0.76
    for k in range(n + 1):
        v = k*h
        ax.plot([v, v], [0, 1], color=LINE, lw=0.25); ax.plot([0, 1], [v, v], color=LINE, lw=0.25)
    m = int(round((hi - lo)/h))
    for k in range(m):
        v = lo + (k + 0.5)*h
        ax.plot([v, v], [lo, hi], color=LINE, lw=0.25); ax.plot([lo, hi], [v, v], color=LINE, lw=0.25)
    ax.add_patch(Rectangle((lo, lo), hi - lo, hi - lo, fc="none", ec="#eb6834", lw=0.8))
    water(ax, np.array([[0, 0], [1, 0], [1, 1], [0, 1]]), eta, (0, 1)); ax.set_title(title)


refined(axes[2], 25, lambda x: 0.49 + 0*x, "(c) refined box")
refined(axes[3], 50, lambda x: 0.5 + 0.05*np.cos(np.pi*x), "(d) standing wave")
# (e) progressive wave
ax = axes[4]
c = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "periodicWaveFlume", "gravityModelScan", "0029")
W = K.Wave(c); M = K.Mesh(c)
L = W.L; th0 = lambda x: W.k*(x - L/2)
water(ax, np.array([[0, M.z1], [L, M.z1], [L, M.z2], [0, M.z2]]), lambda x: W.eta(th0(x)) + M.z1, (0, L))
ax.set_title("(e) progressive wave (periodic)")
ax.set_xlim(-0.02, L + 0.02); ax.set_ylim(M.z1 - 0.03, M.z2 + 0.03)
for ax in axes[:4]:
    ax.set_xlim(-0.08, 1.08); ax.set_ylim(-0.1, 1.1)
for ax in axes:
    ax.set_aspect("equal"); ax.axis("off")
fig.savefig("figures/cases.pdf", bbox_inches="tight"); fig.savefig("figures/cases.png", dpi=150, bbox_inches="tight")
