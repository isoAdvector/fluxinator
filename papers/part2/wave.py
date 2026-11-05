"""Progressive wave (H 0.125 m, T 2 s, d 0.4 m, periodic, one wavelength), interIsoGoam with
the geometric gravity model, Euler, ddtCorr off, one outer corrector unless stated, 10 s.
Per case and convection variant, over 8-10 s: dH (wave height change), r0 and emax against the
exact input wave (hmatch.analyse_input), surface current (wavelength-averaged u - u_exact in
the top full water cell), celerity against the input wave, execution time.
Writes tab_wave/wave.txt and fig_waveCrest/waveCrest.pdf (run by their make.py)."""
import os
import re
import numpy as np
from common import CONV, K, plt, save, write, FULL, COL
from hmatch import analyse_input, celerity

CASES = {
    "190x40, dt 0.005 s": {
        "implicit": "implicit_190x40_dt0.005", "Crank-Nicolson": "crankNicolson_190x40_dt0.005",
        "explicit": "explicit_190x40_dt0.005",
        "explicit + LW": "explicitLW_190x40_dt0.005",
        "explicit + interfaceUpwind": "explicitIU_190x40_dt0.005",
        "explicit + LW + interfaceUpwind": "explicitLWIU_190x40_dt0.005"},
    "190x40, dt 0.01 s": {
        "implicit": "implicit_190x40_dt0.01", "Crank-Nicolson": "crankNicolson_190x40_dt0.01",
        "explicit": "explicit_190x40_dt0.01", "explicit + LW": "explicitLW_190x40_dt0.01",
        "explicit + interfaceUpwind": "explicitIU_190x40_dt0.01",
        "explicit + LW + interfaceUpwind": "explicitLWIU_190x40_dt0.01"},
    "380x80, dt 0.0025 s": {
        "implicit": "implicit_380x80_dt0.0025", "Crank-Nicolson": "crankNicolson_380x80_dt0.0025",
        "explicit": "explicit_380x80_dt0.0025",
        "explicit + LW": "explicitLW_380x80_dt0.0025",
        "explicit + interfaceUpwind": "explicitIU_380x80_dt0.0025",
        "explicit + LW + interfaceUpwind": "explicitLWIU_380x80_dt0.0025"}}


def surface_current(c, t0=8, t1=10):
    M, W = K.Mesh(c), K.Wave(c); nx, nz = M.nx, M.nz
    xcell = (np.arange(nx) + 0.5)*M.dx; zc = (np.arange(nz) + 0.5)*M.dz
    e = []
    for t in K.times(c):
        if not t0 < float(t) <= t1 + 1e-9:
            continue
        A = K.foamio.read_vol(f"{c}/{t}/alpha.water").reshape(nz, nx)
        U = K.foamio.read_vol(f"{c}/{t}/U")[:, 0].reshape(nz, nx)
        h = A.sum(0)*M.dz; sc, _ = K._parab(h, int(np.argmax(h)), 1); xc = sc*M.dx
        for i in range(nx):
            kW = np.where(A[:, i] > 0.99)[0].max()
            e.append(U[kW, i] - W.uv(W.k*(xcell[i] - xc), zc[kW])[0][0])
    return np.mean(e)


def exectime(c):
    log = [f for f in os.listdir(c) if f.startswith("log.interIsoGoam")][0]
    return float(re.findall(r"^ExecutionTime = ([0-9.]+)", open(f"{c}/{log}").read(), re.M)[-1])


if __name__ == "__main__":
    lines = ["Progressive wave, 8-10 s. dH [%]; r0 = u/u_exact in the top full water cell under the crest and",
             "emax [% of exact crest-surface speed], both against the exact input wave; surface current",
             "[m/s]; celerity error against the input wave [%]; execution time [s] (serial).",
             "Sources: run/progressiveWave/convection/<case>", ""]
    for title, runs in CASES.items():
        lines.append(f"-- {title}")
        lines.append(f"{'variant':34s} {'case':32s} {'dH':>6s} {'r0':>5s} {'emax':>6s} {'current':>8s} {'celer.':>7s} {'time':>7s}")
        for v, c in runs.items():
            p = os.path.join(CONV, c)
            if not os.path.exists(p) or not any(f.startswith("log.") and "End" in open(os.path.join(p, f)).read()
                                                for f in os.listdir(p) if f.startswith("log.interIsoGoam")):
                lines.append(f"{v:34s} {c:32s}  (not available)"); continue
            a, cel = analyse_input(p, (8, 10)), celerity(p)
            lines.append(f"{v:34s} {c:32s} {100*a['dH']:+6.1f} {a['r0']:5.2f} {100*a['emax']:6.1f} "
                         f"{surface_current(p):+8.3f} {100*cel['dc']:+7.2f} {exectime(p):7.0f}")
        lines.append("")
    write("wave.txt", lines)

    # crest profiles at t = 10 s
    show = ["implicit", "explicit", "explicit + interfaceUpwind"]
    fig, axes = plt.subplots(1, 2, figsize=(FULL, 3.3), sharey=True)
    for ax, title in zip(axes, ["190x40, dt 0.005 s", "380x80, dt 0.0025 s"]):
        for v in show:
            p = os.path.join(CONV, CASES[title][v])
            if not os.path.exists(os.path.join(p, "10")):
                continue
            M, W0 = K.Mesh(p), K.Wave(p)
            z, un, _ = K.snapshot(p, "10", M, W0, None)["prof"]
            # distinct marker shapes so the runs can be told apart without colour; interfaceUpwind
            # as open circles without a line, so the exact line can be seen through them
            st = {"implicit": dict(marker="s", ls="-"), "explicit": dict(marker="^", ls="-"),
                  "explicit + interfaceUpwind": dict(marker="o", ls="none", mfc="none", mew=1.1)}[v]
            ax.plot(un, z, color=COL[v], ms=4, lw=1.0, label=v, zorder=3, **st)
        # exact profile as a thin black line on top
        ze = np.linspace(M.z1, M.z1 + W0.eta(0.0)[0], 200)
        ax.plot(W0.uv(0.0, ze - M.z1)[0], ze, "-", color="k", lw=0.9,
                label="exact (input wave)", zorder=5)
        ax.set_title(title.replace("dt", r"$\Delta t$").replace("x", r"$\times$", 1) + ", t = 10 s", fontsize=9)
        ax.set_xlabel("horizontal velocity under the crest [m/s]"); ax.set_ylim(-0.42, 0.11)
    axes[0].set_ylabel("z above still water level [m]")
    axes[0].legend(loc="lower right", fontsize=8)
    fig.tight_layout()
    save(fig, "waveCrest")
