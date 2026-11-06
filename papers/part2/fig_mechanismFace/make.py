# Copyright (C) 2026 Johan Roenby
# SPDX-License-Identifier: GPL-3.0-or-later
"""Figure 6 (fig:mechanismFace): how far one explicit convection step can amplify velocities in the
interface cells of a wave crest, for upwind and non-upwind face values. Uses the crest model of
crestFaces.py / crestSchemes.py (linear two-phase wave, 190x40 mesh, dt 0.005 s, exact mass fluxes).
Gain of a cell = sum over k of |dU_i^{n+1}/dU_k^n|: 1 for a weighted average of old velocities.
Part II, Figure fig:mechanismFace. Writes gain.txt and mechanismFace.pdf here (cache
../data/gain_dt0.005.npz). Usage: python3 make.py"""
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

import os
import sys
import numpy as np
sys.argv = [sys.argv[0], "0.005"]
import crestSchemes as S  # noqa: E402
from crestFaces import L  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, "..", "data", "gain_dt0.005.npz")
SCH = ["interfaceUpwind", "limitedLinearV", "linear"]
cells = [tuple(c) for c in np.argwhere(S.iface)]
if os.path.exists(CACHE):
    G = dict(np.load(CACHE))
else:
    G = {}
    for k in SCH:
        g = S.gain(k, cells)
        G[k] = np.array([g[c] for c in cells])
    np.savez(CACHE, **G)
out = np.zeros((S.NZ, S.NX))
np.add.at(out, (S.fP[:, 0], S.fP[:, 1]), np.clip(S.F, 0, None))
np.add.at(out, (S.fN[:, 0], S.fN[:, 1]), np.clip(-S.F, 0, None))
ratio = np.array([out[c]/S.m1[c] for c in cells])
cn = np.array(cells)
xcell = 0.5*(S.XE[cn[:, 1]] + S.XE[cn[:, 1] + 1]) - S.XC
xcell = (xcell + 0.5*L) % L - 0.5*L          # crest at 0, wave travelling to +x
falling = xcell < 0

if __name__ == "__main__":
    from common import plt, save, write, FULL, COL
    from crestFaces import eta, theta
    lines = ["Gain of one explicit convection step in the interface cells of the crest model (crestFaces.py),",
             "190x40, dt 0.005 s: gain = sum_k |dU_i^{n+1}/dU_k^n| (x component); F_out/m^{n+1} = mass leaving",
             "the cell in the step over its mass after the step. Cells where the surface falls: x - x_c < 0.",
             "%d interface cells, max F_out/m^{n+1} %.1f; F_out/m^n max %.2f" % (len(cells), ratio.max(),
             max(out[c]/S.m0[c] for c in cells)),
             "scheme            max gain  max gain (rising surface)  max |gain - (1 + F_out/m^{n+1})|"]
    for k in SCH:
        lines.append("%-16s %9.2f %10.2f %12.2f" % (k, G[k].max(), G[k][~falling].max(),
                                                   np.abs(G[k] - 1 - ratio).max()))
    write("gain.txt", lines)

    STY = {"interfaceUpwind": dict(color=COL["explicit + interfaceUpwind"], marker="o", ms=3.5, label="interfaceUpwind"),
           "limitedLinearV": dict(color=COL["explicit"], marker="s", ms=3.5, label="limitedLinearV"),
           "linear": dict(color=COL["implicit"], marker="^", ms=3.5, label="linear")}
    fig = plt.figure(figsize=(FULL, 2.7))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.6, 1], height_ratios=[1, 2.6], hspace=0.08, wspace=0.28)
    ax0 = fig.add_subplot(gs[0, 0]); ax = fig.add_subplot(gs[1, 0], sharex=ax0); bx = fig.add_subplot(gs[:, 1])
    x = np.linspace(-0.5*L, 0.5*L, 600)
    ax0.plot(x, eta(theta(x + S.XC, 0.0)), color=COL["exact"])
    ax0.axvspan(-0.5*L, 0, color="#dddddd", lw=0)
    ax0.text(-0.5*L + 0.05, 0.035, "surface falls", fontsize=8, va="center")
    ax0.text(0.05, 0.035, "surface rises", fontsize=8, va="center")
    ax0.annotate("", (0.9, -0.03), (0.4, -0.03), arrowprops=dict(arrowstyle="-|>", color=COL["exact"], lw=0.8))
    ax0.set_ylabel(r"$\eta$ [m]"); ax0.set_ylim(-0.08, 0.08); ax0.tick_params(labelbottom=False)
    ax0.set_title("(a) gain along the wave", fontsize=9, loc="left")
    ax.axvspan(-0.5*L, 0, color="#dddddd", lw=0)
    for k in SCH[::-1]:
        st = dict(STY[k]); ax.plot(xcell, G[k], ls="none", mfc="none", mew=0.8, **st)
    ax.set_xlim(-0.5*L, 0.5*L); ax.set_xlabel(r"$x - x_c$ [m]"); ax.set_ylabel("gain")
    ax.set_ylim(0, 15)
    r = np.linspace(0, ratio.max()*1.05, 50)
    bx.plot(r, 1 + r, color="#888888", lw=0.8, ls="--", label=r"$1 + F_\mathrm{out}/m^{n+1}$")
    for k in SCH[::-1]:
        bx.plot(ratio, G[k], ls="none", mfc="none", mew=0.8, **STY[k])
    bx.set_xlabel(r"$F_\mathrm{out}/m^{n+1}$"); bx.set_ylabel("gain"); bx.set_ylim(0, 15); bx.set_xlim(0, 13.5)
    bx.set_title("(b) gain and mass leaving the cell", fontsize=9, loc="left")
    h, l = bx.get_legend_handles_labels()
    bx.legend(h[::-1], l[::-1], fontsize=8, loc="upper left")
    save(fig, "mechanismFace")
