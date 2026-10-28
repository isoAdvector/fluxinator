"""Progressive wave, 380x80, max Co 0.2, t = 30 s (gravityModelScan 0027 interFoam, 0028
interIsoFoam, 0029 interIsoGoam): error of the horizontal velocity in the water relative to the
exact input wave at the crest position of the run; air (alpha < 0.5) not shown.
Writes figures/waveField.pdf. Usage: python3 figWaveField.py"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "qwistWave"))
import kinematics as K  # noqa: E402
plt.rcParams.update({"font.family": "serif", "font.size": 8})
SCAN = os.path.join(HERE, "..", "periodicWaveFlume", "gravityModelScan")
RUNS = [("0027", "(a) interFoam"), ("0028", "(b) interIsoFoam"), ("0029", "(c) interIsoGoam")]
t = "30"
fig, axes = plt.subplots(3, 1, figsize=(7.0, 4.6), sharex=True)
for ax, (case, name) in zip(axes, RUNS):
    c = os.path.join(SCAN, case); M, W = K.Mesh(c), K.Wave(c); nx, nz = M.nx, M.nz
    A = K.foamio.read_vol(f"{c}/{t}/alpha.water").reshape(nz, nx)
    U = K.foamio.read_vol(f"{c}/{t}/U").reshape(nz, nx, 3)
    s = K.snapshot(c, t, M, W, None)
    x = (np.arange(nx) + 0.5)*M.dx; zc = (np.arange(nz) + 0.5)*M.dz
    X, Z = np.meshgrid(x, zc)
    ue = W.uv(W.k*(X - s["xc"]).ravel(), Z.ravel())[0].reshape(nz, nx)
    ax.grid(False)
    im = ax.pcolormesh(np.arange(nx + 1)*M.dx, M.z1 + np.arange(nz + 1)*M.dz,
                       np.ma.masked_where(A < 0.5, U[:, :, 0] - ue), cmap="RdBu_r", vmin=-0.15, vmax=0.15,
                       shading="flat")
    ax.contour(x, M.z1 + zc, A, [0.5], colors="#1f1f1e", linewidths=0.7)
    ax.set_ylim(-0.4, 0.12); ax.set_ylabel("z [m]")
    ax.set_title(name, fontsize=8.5, loc="left")
axes[-1].set_xlabel("x [m]")
fig.colorbar(im, ax=axes, label=r"$u - u_{\mathrm{exact}}$ [m/s]", shrink=0.8, pad=0.02)
fig.savefig("figures/waveField.pdf", bbox_inches="tight"); fig.savefig("figures/waveField.png", dpi=150, bbox_inches="tight")
