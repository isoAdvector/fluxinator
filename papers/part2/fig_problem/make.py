# Copyright (C) 2026 Johan Roenby
# SPDX-License-Identifier: GPL-3.0-or-later
"""The problem: progressive wave (190x40, dt 0.005 s) with the hardcoded implicit convection.
Horizontal velocity in the column nearest (a) the crest and (b) the trough after one, three and
five periods, against the exact input wave. Under the crest the top cells become too fast; under
the trough the top cells also gain forward velocity and the deeper water slows down: a spurious
current near the surface along the whole wave, with a return flow below.
Part II, Figure fig:problem. Output: problem.pdf. Usage: python3 make.py"""
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

import os
import numpy as np
import matplotlib
from matplotlib import patheffects
from common import CONV, K, plt, save, FULL, COL

c = os.path.join(CONV, "implicit_190x40_dt0.005")
M, W0 = K.Mesh(c), K.Wave(c)
# interIsoGoam with implicit convection: its navy (both papers), lighter for earlier times as the
# successive crests of Part I, mixed with white so that overlapping lines stay opaque
NAVY = np.array(matplotlib.colors.to_rgb(COL["implicit"]))
TIMES = [(t, tuple(1 - f*(1 - NAVY)), mk) for t, f, mk in (("2", 0.35, "o"), ("6", 0.65, "s"), ("10", 1.0, "^"))]


def profiles(t):
    """Wet-cell profiles of u in the columns nearest the fitted crest and trough."""
    A = K.foamio.read_vol(f"{c}/{t}/alpha.water").reshape(M.nz, M.nx)
    U = K.foamio.read_vol(f"{c}/{t}/U")[:, 0].reshape(M.nz, M.nx)
    h = M.z1 + A.sum(0)*M.dz
    z = M.z1 + (np.arange(M.nz) + 0.5)*M.dz
    out = {}
    for name, i in (("crest", int(np.argmax(h))), ("trough", int(np.argmin(h)))):
        wet = A[:, i] >= 0.5
        out[name] = (U[wet, i], z[wet])
    return out


fig, axes = plt.subplots(1, 2, figsize=(FULL, 3.2), sharey=True)
P = {t: profiles(t) for t, _, _ in TIMES}
for ax, name, theta, title in [(axes[0], "crest", 0.0, "(a) under the crest"),
                               (axes[1], "trough", np.pi, "(b) under the trough")]:
    for t, colour, mk in TIMES:
        u, z = P[t][name]
        ax.plot(u, z, "-" + mk, color=colour, label=f"interIsoGoam, t = {t} s")
    ze = np.linspace(M.z1, M.z1 + W0.eta(theta)[0], 200)
    # exact wave on top, thick with long dashes and a white outline, so that it stays visible in
    # poor print where it crosses the profiles
    ax.plot(W0.uv(theta, ze - M.z1)[0], ze, color=COL["exact"], lw=1.6, ls=(0, (5, 2)), zorder=5,
            path_effects=[patheffects.withStroke(linewidth=3.0, foreground="white")], label="exact (input wave)")
    ax.set_xlabel("horizontal velocity [m/s]")
    ax.set_title(title, fontsize=9, loc="left")
axes[0].set_ylabel("z above still water level [m]")
axes[0].set_ylim(-0.42, 0.11)
axes[0].legend(fontsize=8, loc="lower right")
fig.tight_layout()
save(fig, "problem")
