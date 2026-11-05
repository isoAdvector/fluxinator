"""Error of the horizontal velocity in the water of the progressive wave at t = 10 s (190x40,
dt 0.005 s): u - u_exact, with the exact input wave at the crest position of the run; air (alpha < 0.5) masked. Implicit convection and explicit convection with
LW + interfaceUpwind. Part II, Figure fig:waveField. Output: waveField.pdf. Usage: python3 make.py"""
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

import os
import numpy as np
from common import CONV, K, plt, save, FULL

RUNS = [("implicit_190x40_dt0.005", "(a) implicit convection"),
        ("explicitIU_190x40_dt0.005", "(b) explicit convection, interfaceUpwind")]
fig, axes = plt.subplots(2, 1, figsize=(FULL, 3.6), sharex=True)
for ax, (c, title) in zip(axes, RUNS):
    c = os.path.join(CONV, c); t = "10"
    M = K.Mesh(c); nx, nz = M.nx, M.nz
    W = K.Wave(c); s = K.snapshot(c, t, M, W, None)
    A = K.foamio.read_vol(f"{c}/{t}/alpha.water").reshape(nz, nx)
    U = K.foamio.read_vol(f"{c}/{t}/U")[:, 0].reshape(nz, nx)
    x = (np.arange(nx) + 0.5)*M.dx; z = (np.arange(nz) + 0.5)*M.dz
    Xc, Zc = np.meshgrid(x, z)
    ue = W.uv(W.k*(Xc - s["xc"]).ravel(), Zc.ravel())[0].reshape(nz, nx)
    e = np.ma.masked_where(A < 0.5, U - ue)
    ax.grid(False)
    im = ax.pcolormesh(np.arange(nx + 1)*M.dx, M.z1 + np.arange(nz + 1)*M.dz, e, cmap="RdBu_r",
                       vmin=-0.08, vmax=0.08, shading="flat")
    ax.contour(x, M.z1 + z, A, [0.5], colors="#1f1f1e", linewidths=0.7)
    ax.set_ylim(-0.4, 0.15); ax.set_ylabel("z [m]")
    ax.set_title(title, fontsize=9, loc="left")
axes[1].set_xlabel("x [m]")
fig.colorbar(im, ax=axes, label=r"$u - u_{\mathrm{exact}}$ [m/s]", shrink=0.9)
save(fig, "waveField")
