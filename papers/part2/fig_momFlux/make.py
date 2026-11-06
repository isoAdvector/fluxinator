# Copyright (C) 2026 Johan Roenby
# SPDX-License-Identifier: GPL-3.0-or-later
"""Wave-averaged vertical flux of horizontal momentum below the wave trough, 190x40 progressive
wave, 2 < t <= 10 s: solver value rho_f phi_f U_f,x / A_f on the horizontal faces, with U_f the
face velocity the convection term used (theta-weighted between U^n and U^{n+1}; function object
convectiveMomentumFlux), against the prediction rho_w <w(t) u(t + theta dt)> of the exact wave at
fixed height (zero for theta = 0). Runs: run/progressiveWave/momentumFlux.
Part II, Figure fig:momFlux. Output: momFlux.pdf, momFlux.txt. Usage: python3 make.py"""
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

import os
import re
import numpy as np
from common import WAVE, K, plt, save, write, HALF, COL

MF = os.path.join(WAVE, "momentumFlux")

RHO_W = 1000.0
RUNS = [("implicit_190x40_dt0.005", 1.0, "implicit, $\\Delta t$ = 0.005 s", COL["implicit"], "-"),
        ("implicit_190x40_dt0.0025", 1.0, "implicit, $\\Delta t$ = 0.0025 s", COL["implicit"], "--"),
        ("crankNicolson_190x40_dt0.005", 0.5, "Crank-Nicolson, $\\Delta t$ = 0.005 s", COL["Crank-Nicolson"], "-"),
        ("explicit_190x40_dt0.005", 0.0, "explicit, $\\Delta t$ = 0.005 s", COL["explicit"], "-")]


def profile(case):
    M = K.Mesh(case); nx, nz = M.nx, M.nz
    pm = f"{case}/constant/polyMesh"
    own = K.foamio.read_labels(f"{pm}/owner"); nei = K.foamio.read_labels(f"{pm}/neighbour")
    nI = nei.size; own = own[:nI]
    hf = np.where(nei - own == nx)[0]; ko = own[hf]//nx
    dt = float(re.search(r"^deltaT\s+([0-9.eE+-]+);", open(f"{case}/system/controlDict").read(), re.M).group(1))
    S, n = np.zeros(nz), 0
    for t in K.times(case):
        if not 2 < float(t) <= 10 + 1e-9:
            continue
        f = K.foamio.read_vol(f"{case}/{t}/convectiveMomentumFlux")[:nI, 0][hf]
        np.add.at(S, ko, f/(M.dx*0.2)); n += 1
    z = M.z1 + (np.arange(nz) + 1)*M.dz
    return z, S/(n*nx), dt


def predicted(case, z, lag):
    W = K.Wave(case); M = K.Mesh(case)
    th = np.linspace(0, 2*np.pi, 721)[:-1]
    out = []
    for zz in z:
        u = W.uv(th - W.omega*lag, np.full_like(th, zz - M.z1))[0]
        w = W.uv(th, np.full_like(th, zz - M.z1))[1]
        out.append(RHO_W*np.mean(u*w))
    return np.array(out)


lines = ["Wave-averaged upward flux of x-momentum [N/m^2] below the trough (trough at -0.043 m), 190x40, 2-10 s.",
         "solver: rho_f phi_f U_f,x/A_f with the face velocity used; predicted: rho_w <w(t) u(t + theta dt)>.",
         "Sources: run/progressiveWave/momentumFlux/<case> (convectiveMomentumFlux)", ""]
fig, ax = plt.subplots(figsize=(1.6*HALF, 2.9))
for case, theta, lab, c, ls in RUNS:
    z, S, dt = profile(os.path.join(MF, case))
    sel = (z > -0.2) & (z < -0.05)
    P = predicted(os.path.join(MF, case), z[sel], theta*dt)
    ax.plot(S[sel], z[sel], ls, color=c, marker="o", label=lab)
    if theta > 0:
        ax.plot(P, z[sel], ls, color=c, lw=0.7, alpha=0.7)
    lines.append(f"{case} (theta {theta}, dt {dt}):")
    lines += [f"   z {zz:+.3f}  solver {s:+.4f}  predicted {p:+.4f}" for zz, s, p in zip(z[sel], S[sel], P)]
ax.axvline(0, color="#1f1f1e", lw=0.6)
ax.plot([], [], "-", color="#888888", lw=0.7, label="thin lines: predicted error")
ax.set_xlabel(r"Momentum flux $\rho \langle w(t)\, u(t+\theta\Delta t)\rangle$ [N/m$^2$]"); ax.set_ylabel("z above still water level [m]")
h, l = ax.get_legend_handles_labels()
ax.legend(h, l, loc="upper left", bbox_to_anchor=(1.03, 1.0), ncol=1, fontsize=8, borderaxespad=0)
fig.tight_layout()
save(fig, "momFlux")
write("momFlux.txt", lines)
