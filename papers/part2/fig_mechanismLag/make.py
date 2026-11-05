"""Schematic of the two mechanisms. (a) Time lag: at a fixed point below the trough of a
progressive wave, u and w are 90 degrees out of phase, so <w u> = 0; the convection term pairs the mass
flux at t with the velocity at t + theta*dt, and <w(t) u(t + dt)> > 0 (dt exaggerated, T/10).
(b) Face value at an emptying interface cell: water leaves cell P into the air cell N; with the
upwind face value the momentum leaves with U_P and U_P is unchanged; with a central value it
leaves with (U_P + U_N)/2, and the momentum difference is divided by the small mass left in P.
Part II, Figure fig:mechanismLag. Output: mechanismLag.pdf. Usage: python3 make.py"""
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

import numpy as np
from matplotlib.patches import Rectangle, FancyArrowPatch
from common import plt, save, FULL, COL

# Panels side by side, both with the size of (a) when the panels were stacked
# (0.55*FULL x 3.3 in figure): 2.685 x 1.355 in
W, H, x0, gap, y0 = 6.9, 2.15, 0.62, 0.78, 0.45
fig = plt.figure(figsize=(W, H))
axes = [fig.add_axes([x0/W, y0/H, 2.685/W, 1.355/H]),
        fig.add_axes([(x0 + 2.685 + gap)/W, y0/H, 2.685/W, 1.355/H])]
t = np.linspace(0, 1, 400); lag = 0.1
w, u, ul = -np.sin(2*np.pi*t), np.cos(2*np.pi*t), np.cos(2*np.pi*(t + lag))
ax = axes[0]
ax.plot(t, w, color=COL["exact"], label="$w(t)$")
ax.plot(t, u, "--", color=COL["implicit"], label="$u(t)$")
ax.plot(t, ul, color=COL["implicit"], label=r"$u(t+\Delta t)$")
ax.set_xlabel("t / T"); ax.set_ylabel("velocity / amplitude")
ax.set_title("(a) velocities", fontsize=9, loc="left")
ax.legend(fontsize=8, ncol=3, loc="lower center", columnspacing=1.2, handlelength=1.6); ax.set_ylim(-2.1, 1.2)
ax = axes[1]
ax.plot(t, w*u, "--", color="#888888", label=r"$w(t)\,u(t)$, mean 0")
ax.plot(t, w*ul, color=COL["implicit"], label=r"$w(t)\,u(t+\Delta t)$, mean $\neq 0$")
ax.axhline(np.mean(w*ul), color=COL["implicit"], lw=0.7, ls=":")
ax.set_xlabel("t / T"); ax.set_ylabel("product"); ax.legend(fontsize=8, loc="lower left")
ax.set_title("(b) momentum flux", fontsize=9, loc="left")
ax.set_ylim(-1.55, 1.1)
save(fig, "mechanismLag")

# mechanismFace (Figure 6) is now made by figGain.py from the crest model
