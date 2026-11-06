# Copyright (C) 2026 Johan Roenby
# SPDX-License-Identifier: GPL-3.0-or-later
"""Full solver (interIsoGoam) on the periodic disc case (run/periodicDisc) (tutorial discInConstantFlowCyclicBCs,
50x50, g = 0, U0 = (0.1, 0, 0.05) m/s, dt 0.005 s, 1 s, ddtCorr off):
(1) uniform acceleration a = (0.5, 0, 0) m/s^2 applied as the body force rho*a with rho^{n+1}
(new) or rho^n (old); (2) uniform flow without force, all convection schemes, density ratio
1e3 (1 s) and 1e4 (3 s). Exact: U = U0 + a t, uniform p_rgh. Writes tab_fullSolver/fullSolver.txt and
fig_fullSolver/fullSolver.pdf (|U - U0 - a t| at t = 1 s), run by their make.py."""
import os
import numpy as np
from common import DISC, K, plt, save, write, FULL

CT = DISC
U0 = np.array([0.1, 0, 0.05])


def errors(case, a):
    m = np.zeros(4)
    for t in K.times(case):
        if float(t) == 0 or not os.path.exists(f"{case}/{t}/U"):
            continue
        A = K.foamio.read_vol(f"{case}/{t}/alpha.water"); U = K.foamio.read_vol(f"{case}/{t}/U")
        p = K.foamio.read_vol(f"{case}/{t}/p_rgh")
        e = np.linalg.norm(U - (U0 + a*float(t)), axis=1)
        w, i = A > 0.99, (A >= 0.01) & (A <= 0.99)
        m = np.maximum(m, [e[w].max(), e[i].max() if i.any() else 0, e[~(w | i)].max(), p.max() - p.min()])
    return m


acc = np.array([0.5, 0, 0])
rows = [("implicit_rhoNew", "implicit, limitedLinearV, rho^{n+1}"), ("implicit_rhoOld", "implicit, limitedLinearV, rho^n"),
        ("implicitIU_rhoNew", "implicit, interfaceUpwind, rho^{n+1}"), ("explicit_rhoNew", "explicit, limitedLinearV, rho^{n+1}"),
        ("explicit_rhoOld", "explicit, limitedLinearV, rho^n"), ("explicitIU_rhoNew", "explicit, interfaceUpwind, rho^{n+1}"),
        ("explicitLW_rhoNew", "explicit, LW, rho^{n+1}"), ("explicitLWIU_rhoNew", "explicit, LW + interfaceUpwind, rho^{n+1}")]
lines = ["Full solver, uniform acceleration a = (0.5 0 0) m/s^2, body force rho^{n+1} a (new) or rho^n a (old), system/fvOptions.",
         "max over t <= 1 s of |U - U0 - a t| [m/s] in water (alpha > 0.99), interface, air; p_rgh range [Pa].",
         "Sources: run/periodicDisc/<case> (runCase.sh, baseCase)", "",
         f"{'case':20s} {'variant':42s} {'water':>9s} {'interf.':>9s} {'air':>9s} {'p_rgh':>9s}"]
for c, lab in rows:
    m = errors(os.path.join(CT, c), acc)
    lines.append(f"{c:20s} {lab:42s} " + " ".join(f"{x:9.1e}" for x in m))
lines += ["", "Uniform flow without force (Liu et al. 2024, Fig. 10): density ratio 1e3 (noForce_*, 1 s) and 1e4 (noForceRatio1e4_*, 3 s)",
          f"{'case':32s} {'water':>9s} {'interf.':>9s} {'air':>9s} {'p_rgh':>9s}"]
for c in ["noForce_implicitUpwind", "noForce_implicit", "noForce_implicitLinear", "noForce_implicit_ddtCorrOn", "noForce_explicit", "noForce_explicitIU",
          "noForce_explicitLinear", "noForceRatio1e4_implicitUpwind", "noForceRatio1e4_implicit", "noForceRatio1e4_implicitLinear", "noForceRatio1e4_implicitQUICK", "noForceRatio1e4_implicitCubic", "noForceRatio1e4_explicitIU"]:
    m = errors(os.path.join(CT, c), np.zeros(3))
    lines.append(f"{c:32s} " + " ".join(f"{x:9.1e}" for x in m))
write("fullSolver.txt", lines)

n = 50
x = (np.arange(n) + 0.5)/n
fig, axes = plt.subplots(1, 2, figsize=(FULL*0.8, 3.0), sharey=True)
for ax, (c, lab) in zip(axes, [("implicit_rhoNew", "(a) implicit convection"), ("explicit_rhoNew", "(b) explicit convection")]):
    case = os.path.join(CT, c)
    e = np.linalg.norm(K.foamio.read_vol(f"{case}/1/U") - (U0 + acc), axis=1).reshape(n, n)
    A = K.foamio.read_vol(f"{case}/1/alpha.water").reshape(n, n)
    ax.grid(False)
    # each panel its own linear colour scale, with the colour bar inside the panel (as in Part I)
    im = ax.pcolormesh(np.linspace(0, 1, n + 1), np.linspace(0, 1, n + 1), e,
                       vmin=0, vmax=e.max(), cmap="magma_r", shading="flat")
    ax.contour(x, x, A, [0.5], colors="#2a78d6", linewidths=0.8)
    cax = ax.inset_axes([0.08, 0.15, 0.52, 0.045], zorder=4)
    cax.grid(False)
    cb = fig.colorbar(im, cax=cax, orientation="horizontal", ticks=[0, e.max()])
    cb.ax.set_xticklabels(["0", f"{e.max():.1e}"], fontsize=8)
    cb.ax.tick_params(length=2, pad=1)
    cax.set_title(r"$|\mathbf{U} - \mathbf{U}_0 - \mathbf{a}t|$ [m/s]", fontsize=8, pad=2)
    # direction of the flow, U0 + a t at t = 1 s
    d = (U0 + acc)[[0, 2]]/np.linalg.norm((U0 + acc)[[0, 2]])
    ax.annotate("", xy=(0.08 + 0.3*d[0], 0.86 + 0.3*d[1]), xytext=(0.08, 0.86), xycoords="axes fraction",
                arrowprops=dict(arrowstyle="-|>", color="k", lw=1.2, mutation_scale=10), zorder=5)
    ax.text(0.08, 0.89, "flow", transform=ax.transAxes, fontsize=8, va="bottom", zorder=5)
    ax.set_aspect("equal"); ax.set_title(lab, fontsize=9, loc="left")
    ax.set_xticks([]); ax.set_yticks([])
fig.tight_layout(w_pad=1.0)
# white box fitted around each colour bar with its title and tick labels
fig.canvas.draw()
for ax in axes:
    cax = ax.child_axes[0]
    bb = cax.get_tightbbox(fig.canvas.get_renderer()).transformed(ax.transAxes.inverted())
    pad = 0.025
    ax.add_patch(plt.Rectangle((bb.x0 - pad, bb.y0 - pad), bb.width + 2*pad, bb.height + 2*pad,
                               transform=ax.transAxes, facecolor="white", edgecolor="0.6", lw=0.5,
                               alpha=0.9, zorder=3))
save(fig, "fullSolverUerr")
