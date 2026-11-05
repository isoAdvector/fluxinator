"""The time lag between mass flux and convected velocity: sketch and evidence.
  lagTGV    (a) translating Taylor-Green vortex u = U0 + sin(x - U0 t) cos z, w = -cos(x - U0 t) sin z:
              predicted spurious mean flow <u> - U0 = -theta dt U0 t cos(2z)/2
            (b) toy runs (run/taylorGreenVortex): cos 2z amplitude of <u> - U0 divided by dt U0/2 against t,
              with the predictions -theta t (theta 1 implicit, 1/2 Crank-Nicolson, 0 explicit)
  lagWave   (a) progressive wave, linear theory: lagged vertical flux of horizontal momentum
              F = theta dt <w du/dt> ~ sinh 2kZ, and the mean velocity change it drives: a deficit
              ~ -cosh 2kZ below, the momentum deposited in the top cell (schematic scale)
            (b) progressive wave of Part 1: measured mean current at depth (ulow, data/browser_runs.json,
              Courant-limited, max Co <= 0.2) against the parameter-free prediction
              -theta (30 s/steps) t A^2 omega k <cosh 2kZ>/2 / u_crest, t = 29 s
Writes fig_lagTGV/lagTGV.pdf, fig_lagWave/lagWave.pdf and data/lag.txt (run by their make.py)."""
import json
import os
import re
import sys
import numpy as np
from matplotlib.patches import FancyArrowPatch
from common import plt, STYLE, COL, DATA, FULL, RUN, WAVE, TGV, save, write

sys.path.insert(0, WAVE)
import foamio  # noqa: E402
from analyseScan import U_CREST  # noqa: E402

plt.rcdefaults()
plt.rcParams.update(STYLE)
IMPL, CN, EXPL = COL["implicit"], COL["Crank-Nicolson"], COL["explicit + LW + interfaceUpwind"]
TOY = TGV
ARROW = dict(arrowstyle="-|>", mutation_scale=9, lw=1.0)


def sketch_wave(a):
    # wave column: flux ~ sinh 2kZ (kd ~ 1), deficit ~ -cosh 2kZ, deposit in the top cell
    kd = 1.0
    Z = np.linspace(0, 1, 200)
    F = np.sinh(2*kd*Z)/np.sinh(2*kd)
    du = -np.cosh(2*kd*Z)/np.cosh(2*kd)
    a.plot(F, Z, color=IMPL, lw=1.6)
    for z in (0.25, 0.5, 0.75):
        f = np.sinh(2*kd*z)/np.sinh(2*kd)
        a.add_patch(FancyArrowPatch((f*0.5, z - 0.07), (f*0.5, z + 0.07), color=IMPL, **ARROW))
    m = Z < 0.93
    a.plot(1.55 + 0.45*du[m], Z[m], color="#1f1f1e", lw=1.6)
    a.plot([1.55 + 0.45*du[m][-1], 1.55 + 0.9, 1.55 + 0.9], [Z[m][-1], Z[m][-1], 1.0], color="#1f1f1e", lw=1.6)
    a.axvline(1.55, color="0.5", lw=0.7, ls="--")
    a.axhline(1.0, color="#2a78d6", lw=2.0)
    a.text(0.02, 1.04, "surface: the flux cannot leave", fontsize=8, color="#2a78d6")
    a.text(0.45, -0.05, "lagged upward flux of\nhorizontal momentum\n" r"$\theta\,\Delta t\,\rho\,\langle w\,\partial_t u\rangle$",
           fontsize=8, color=IMPL, va="top", ha="center")
    a.text(2.0, -0.05, "resulting change of\nthe mean velocity", fontsize=8, va="top", ha="center")
    a.text(1.95, 0.84, "surface\ncurrent", fontsize=8, ha="left", va="top")
    a.text(1.02, 0.42, "deficit\nbelow", fontsize=8, ha="right")
    a.set(xlim=(-0.05, 2.65), ylim=(0, 1.13), xticks=[], yticks=[0, 1], yticklabels=["Bed", "SWL"])
    a.grid(False)
    a.set_title("(a) lagged flux and its effect", fontsize=9, loc="left")


def sketch_tgv(b):
    # translating vortex
    x = np.linspace(0, 2*np.pi, 120); z = np.linspace(0, 2*np.pi, 120)
    X, Zz = np.meshgrid(x, z)
    b.contour(X, Zz, np.sin(X)*np.sin(Zz), levels=[-0.75, -0.4, -0.1, 0.1, 0.4, 0.75], colors="0.55",
              linewidths=0.7)
    zz = np.linspace(0, 2*np.pi, 13)[:-1] + np.pi/12
    for q in zz:
        v = -np.cos(2*q)
        b.add_patch(FancyArrowPatch((7.3, q), (7.3 + 0.9*v, q), color=IMPL, clip_on=False, **ARROW))
    zp = np.linspace(0, 2*np.pi, 200)
    b.plot(7.3 - 0.9*np.cos(2*zp), zp, color=IMPL, lw=1.2, clip_on=False)
    b.axvline(7.3, color="0.5", lw=0.7, ls="--", clip_on=False)
    b.text(7.3, -0.75, r"$\langle u\rangle - U_0 \propto -\theta\,\Delta t\,U_0\,t\,\cos 2z$",
           fontsize=8, ha="center", color=IMPL)
    b.set(xlim=(0, 2*np.pi), ylim=(0, 2*np.pi), xticks=[], yticks=[], aspect="equal")
    b.grid(False)
    b.set_title(r"(a) vortex moving with $U_0$ $\rightarrow$", fontsize=9, loc="left")


def amp(case):
    n = 64
    zc = (np.arange(n) + 0.5)*2*np.pi/n
    U0 = float(open(os.path.join(case, "U0")).read())
    ts = sorted((d for d in os.listdir(case) if re.fullmatch(r"\d+(\.\d+)?", d)), key=float)
    t, a = [], []
    for d in ts:
        U = foamio.read_vol(os.path.join(case, d, "U")).reshape(n, n, 3)
        t.append(float(d)); a.append(2*np.mean((U[:, :, 0].mean(1) - U0)*np.cos(2*zc)))
    dt = float(re.search(r"^deltaT\s+(\S+);", open(os.path.join(case, "system/controlDict")).read(), re.M).group(1))
    return np.array(t), np.array(a), dt


def toy(a, rows):
    rows.append("Toy: final amplitude / (dt/2), prediction -theta t_end")
    runs = [("implicit_dt0.02", IMPL, "o", r"implicit, $\Delta t$ 0.02"),
            ("implicit_dt0.01", IMPL, "s", r"implicit, $\Delta t$ 0.01"),
            ("implicit_dt0.005", IMPL, "^", r"implicit, $\Delta t$ 0.005"),
            ("crankNicolson_dt0.01", CN, "D", "Crank–Nicolson"),
            ("explicitLW_dt0.01", CN, "v", "explicit + Lax–Wendroff"),
            ("explicit_dt0.01", EXPL, "o", "explicit"),
            ("implicitNoTranslation_dt0.01", "0.45", "x", r"implicit, $U_0 = 0$")]
    for case, col, mk, lab in runs:
        t, y, dt = amp(os.path.join(TOY, case))
        rows.append(f"{case:18s} a(10)/(dt/2) {y[-1]/(0.5*dt):+7.3f}")
        if "LW" in case:
            continue  # in lag.txt only; LW is introduced after this figure (Section 6)
        a.plot(t, y/(0.5*dt), color=col, marker=mk, ms=3.5, lw=0.9, mfc="white" if mk not in "x" else col,
               label=lab)
    tt = np.array([0, 10])
    for th, col in ((1, IMPL), (0.5, CN), (0, EXPL)):
        a.plot(tt, -th*tt, "--", color=col, lw=1.0)
    a.text(10.2, -10, r"$\theta=1$", fontsize=8, va="center", color=IMPL)
    a.text(10.2, -5, r"$\theta=\frac{1}{2}$", fontsize=8, va="center", color=CN)
    a.text(10.2, 0, r"$\theta=0$", fontsize=8, va="center", color=EXPL)
    a.set(xlabel="t", ylabel=r"$\cos 2z$ amplitude of $\langle u\rangle-U_0$, $/\,\frac{1}{2} U_0\Delta t$",
          xlim=(0, 12), ylim=(-17, 0.8))
    a.legend(fontsize=8, loc="lower left", ncol=1)
    a.set_title("(b) markers: measured; dashed: predicted", fontsize=9, loc="left")


def wave(b, rows):
    # predicted against measured mean current at depth
    H, T, d, L = 0.125, 2.0, 0.4, 3.7873860376294
    am, w, k = H/2, 2*np.pi/T, 2*np.pi/L
    A = am*w/np.sinh(k*d)
    Z = np.linspace(0, 0.5*(d + 0.0708), 200)
    cm = np.mean(np.cosh(2*k*Z))
    R = [r for r in json.load(open(os.path.join(DATA, "browser_runs.json")))
         if r["scan"] == "maxCo" and r["ulow"] is not None and r["par"] <= 0.2 and r["solver"] != "centred"]
    sty = {"interFoam": ("#2a78d6", "s"), "interIsoFoam": ("#eb6834", "^"), "interIsoGoam": ("#1baf7a", "o"),
           "final": ("#1f1f1e", "D")}
    P, M = [], []
    for s, (col, mk) in sty.items():
        th = 0 if s == "final" else 1
        x = [-th*(30.0/r["steps"])*29.0*0.5*A**2*w*k*cm/U_CREST for r in R if r["solver"] == s]
        y = [r["ulow"] for r in R if r["solver"] == s]
        b.plot(x, y, mk, color=col, ms=4.5, mfc="white" if s != "final" else col,
               label=s)
        if s != "final":
            P += x; M += y
    lim = [-0.24, 0.04]
    b.plot(lim, lim, color="0.45", lw=0.8, ls="--")
    b.plot(lim, [2*v for v in lim], color="0.75", lw=0.8, ls=":")
    b.text(-0.088, -0.205, "2 : 1", fontsize=8, color="0.5")
    b.text(-0.17, -0.152, "1 : 1", fontsize=8, color="0.45")
    b.set(xlim=lim, ylim=lim, xlabel=r"predicted mean current at depth $/\,u_\mathrm{crest}$",
          ylabel=r"measured $/\,u_\mathrm{crest}$", aspect="equal")
    b.legend(fontsize=8, loc="upper left")
    b.set_title(r"(b) wave, max Co $\leq$ 0.2, all meshes", fontsize=9, loc="left")
    P, M = np.array(P), np.array(M)
    rows += ["", f"Wave, Part 1 solvers, max Co <= 0.2: measured/predicted median {np.median(M/P):.2f}, "
             f"correlation {np.corrcoef(P, M)[0, 1]:.3f}, n {len(P)}"]


if __name__ == "__main__":
    rows = []
    fig, (a, b) = plt.subplots(1, 2, figsize=(FULL, 2.9), gridspec_kw=dict(width_ratios=(0.8, 1.2)))
    sketch_tgv(a)
    toy(b, rows)
    fig.tight_layout(w_pad=3)
    save(fig, "lagTGV")
    fig, (a, b) = plt.subplots(1, 2, figsize=(FULL, 2.9))
    sketch_wave(a)
    wave(b, rows)
    fig.tight_layout(w_pad=3)
    save(fig, "lagWave")
    write("lag.txt", rows)
