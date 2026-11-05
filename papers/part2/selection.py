"""Selected Courant-limited results for the paper, from data/browser.json (browserData.py).
Selection rule, fixed before looking at the figures: the 380x80 mesh (12 cells per wave height)
at max Co 0.1; the wave height for all four Courant numbers of the same mesh; the cost at max Co 0.1.
  selCase:   380x80, max Co 0.1: (a) wave height per wave, (b) largest speed in the air,
             (c) crest-column velocity profile at 30 s, (d) surface near the crest at 30 s
  selHeight: 380x80, wave height per wave at max Co 0.5, 0.2, 0.1, 0.05
  selCost:   max Co 0.1, 95x20 to 380x80: clock time, time steps and time per step relative to
             interFoam (ddtCorr off). Time steps of the full 30 s runs; time per step from
             run/progressiveWave/data/timing_maxCo0.1.csv (three repeats, no wave gauge); clock
             time = time steps x time per step; whiskers: range of the repeats
Writes fig_sel{Case,Height,Cost}/*.pdf and data/paperSelection.txt (run by their make.py)."""
import csv
import json
import os
import numpy as np
from common import plt, STYLE, DATA, FULL, WAVE, save, write

plt.rcdefaults()
plt.rcParams.update(STYLE)
SOLVERS = ["interFoam", "interIsoFoam", "interIsoGoam", "explicit + LW + interfaceUpwind"]
LABEL = {"interFoam": "interFoam", "interIsoFoam": "interIsoFoam", "interIsoGoam": "interIsoGoam",
         "explicit + LW + interfaceUpwind": "final"}
COLOUR = {"interFoam": "#2a78d6", "interIsoFoam": "#eb6834", "interIsoGoam": "#1baf7a",
          "explicit + LW + interfaceUpwind": "#1f1f1e"}
LW = {s: (1.8 if s.startswith("explicit") else 1.1) for s in SOLVERS}
LS = {s: ((0, (4, 2)) if s.startswith("centred") else "-") for s in SOLVERS}
MESH, CO = 380, 0.1

B = json.load(open(os.path.join(DATA, "browser.json")))
C = {(c["scan"], c["nx"], c["par"], c["solver"]): c for c in B["cases"]}
E = B["exact"]


def get(nx, co, s):
    return C.get(("maxCo", nx, co, s))


def legend(fig, ncol=5, y=1.0):
    h = [plt.Line2D([], [], color=COLOUR[s], lw=LW[s], ls=LS[s], label=LABEL[s]) for s in SOLVERS
         if any(k[3] == s for k in C)]
    h.append(plt.Line2D([], [], color="0.45", lw=1.0, ls="--", label="exact"))
    fig.legend(handles=h, loc="upper center", ncol=ncol, bbox_to_anchor=(0.5, y))


def fig_case():
    fig, axs = plt.subplots(2, 2, figsize=(FULL, 5.2))
    a, b, c, d = axs.flat
    # zoom on the top of the crest, inside (d) at the centre bottom
    dz = d.inset_axes([0.3, 0.04, 0.4, 0.4])
    for s in SOLVERS:
        r = get(MESH, CO, s)
        if r is None:
            continue
        kw = dict(color=COLOUR[s], lw=LW[s], ls=LS[s])
        a.plot(r["tH"], r["dH"], marker="o", ms=2.5, **kw)
        b.plot(r["tair"], r["Uair"], **kw)
        sn = r["snaps"]["30"]
        c.plot(sn["u"], sn["z"], marker="o", ms=2.5, **kw)
        x, h = np.array(sn["x"]), np.array(sn["h"])
        m = np.abs(x) <= 0.12
        d.plot(x[m], h[m], marker="o", ms=2.5, **kw)
        dz.plot(x[m], h[m], marker="o", ms=2.5, **kw)
    ex = np.array(E["x"]); m = np.abs(ex) <= 0.12
    c.plot(E["u"], E["z"], "--", color="0.45", lw=1.0)
    d.plot(ex[m], np.array(E["h"])[m], "--", color="0.45", lw=1.0)
    dz.plot(ex[m], np.array(E["h"])[m], "--", color="0.45", lw=1.0)
    dz.set(xlim=(-0.015, 0.015), ylim=(0.63, 0.665))
    dz.tick_params(labelsize=8, length=2, pad=1, labelbottom=False)
    d.indicate_inset_zoom(dz, edgecolor="0.4")
    a.axhline(0, color="0.45", lw=0.7, ls="--")
    b.axhline(1, color="0.45", lw=0.7, ls="--")
    a.set(xlabel="t [s]", ylabel=r"wave height error [%]", xlim=(0, 30))
    b.set(xlabel="t [s]", ylabel=r"max $|\mathbf{U}|$ in air / $u_\mathrm{crest}$", xlim=(0, 30))
    c.set(xlabel=r"$u_x / u_\mathrm{crest}$", ylabel="height above bed / crest elevation")
    d.set(xlabel="x / wavelength (crest at 0)", ylabel=r"elevation / $H$")
    for ax, t in zip(axs.flat, "abcd"):
        ax.text(0.02, 0.96, f"({t})", transform=ax.transAxes, va="top")
    legend(fig)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    save(fig, "selCase")


def fig_height():
    fig, axs = plt.subplots(1, 4, figsize=(FULL, 2.4), sharex=True)
    for ax, co in zip(axs, (0.5, 0.2, 0.1, 0.05)):
        for s in SOLVERS:
            r = get(MESH, co, s)
            if r is not None:
                ax.plot(r["tH"], r["dH"], color=COLOUR[s], lw=LW[s], ls=LS[s])
        ax.axhline(0, color="0.45", lw=0.7, ls="--")
        ax.set_title(f"max Co {co}", fontsize=9)
        ax.set_xlabel("t [s]")
        ax.set_xlim(0, 30)
    axs[0].set_ylabel("wave height error [%]")
    legend(fig, ncol=5)
    fig.tight_layout(rect=(0, 0, 1, 0.86))
    save(fig, "selHeight")


def fig_cost(co=CO, ref="interFoam"):
    meshes = [95, 190, 380]
    sol = SOLVERS
    T = {}
    for r in csv.DictReader(open(os.path.join(WAVE, "data", "timing_maxCo0.1.csv"))):
        s = SOLVERS[3] if r["solver"] == "final" else r["solver"]
        T.setdefault((int(r["nx"]), s), []).append(float(r["msPerStep"]))
    rows = [f"Relative cost at max Co {co}, normalised by {ref} (clock time, steps, time per step)"]
    fig, axs = plt.subplots(1, 3, figsize=(FULL, 2.5))
    w = 0.8/len(sol)
    for k, (ax, q, lab) in enumerate(zip(axs, ("clock", "steps", "perStep"),
                                         ("clock time", "time steps", "time per step"))):
        for i, s in enumerate(sol):
            vals, lo, hi = [], [], []
            for nx in meshes:
                ms, ms0 = np.array(T[nx, s]), np.mean(T[nx, ref])
                steps = get(nx, co, s)["steps"]/get(nx, co, ref)["steps"]
                x = dict(clock=steps*ms/ms0, steps=np.array([steps]), perStep=ms/ms0)[q]
                vals.append(x.mean())
                lo.append(x.mean() - x.min())
                hi.append(x.max() - x.mean())
                if k == 0:
                    rows.append(f"{nx:4d} {LABEL[s]:18s} clock {x.mean():5.2f} steps {steps:5.2f} "
                                f"perStep {ms.mean()/ms0:5.2f}")
            ax.bar(np.arange(len(meshes)) + (i - (len(sol) - 1)/2)*w, vals, w*0.9,
                   color=COLOUR[s], edgecolor="white", lw=0.5,
                   yerr=None if q == "steps" else [lo, hi],
                   error_kw=dict(lw=0.7, capsize=1.5, ecolor="0.2"))
        ax.axhline(1, color="0.3", lw=0.7)
        ax.set_xticks(range(len(meshes)))
        ax.set_xticklabels([f"{nx}x{nx*4//19}" for nx in meshes], fontsize=8, rotation=20)
        ax.set_title(f"{lab}, relative", fontsize=9)
        ax.grid(axis="x", visible=False)
    h = [plt.Rectangle((0, 0), 1, 1, color=COLOUR[s], label=LABEL[s]) for s in sol]
    fig.legend(handles=h, loc="upper center", ncol=len(sol), bbox_to_anchor=(0.5, 1.0))
    axs[0].set_ylabel(f"relative to {ref}")
    fig.tight_layout(rect=(0, 0, 1, 0.86))
    save(fig, "selCost")
    return rows


if __name__ == "__main__":
    fig_case()
    fig_height()
    rows = fig_cost()
    rows += ["", f"{MESH}x{MESH*4//19}, max Co {CO}: dHfit, dc, Uair median, emax, ulow"]
    for s in SOLVERS:
        r = get(MESH, CO, s)
        if r:
            rows.append(f"{LABEL[s]:18s} dH {r['dHfit']:+6.2f}% dc {r['dc']:+6.3f}% Uair {r['UairMed']:.3f} "
                        f"emax {r['emax']:5.2f}% ulow {r['ulow']:+.4f}")
    write("paperSelection.txt", rows)
