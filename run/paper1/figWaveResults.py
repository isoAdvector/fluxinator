"""Progressive-wave results of Part I as figures, from the metrics of all runs
(run/periodicWaveFlume/paperData/runs_fixedDt.csv and runs_maxCo.csv, written by analyseScan.py).
figures/waveResults.pdf: for each mesh, (top) wave height error after 15 periods (dHfit) and
(bottom) median maximum velocity over t > 2 s relative to the crest velocity (Umed), against the
time-averaged maximum Courant number; filled: fixed time step, open: Courant limit.
figures/waveCost.pdf: runs under a Courant limit, change of interIsoGoam relative to interFoam
(solid) and interIsoFoam (dashed) in (a) the number of time steps and (b) the wall-clock time,
against the Courant limit, one line per mesh.
Usage (in run/paper1): python3 figWaveResults.py"""
import csv
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
PD = os.path.join(HERE, "..", "periodicWaveFlume", "paperData")
plt.rcParams.update({"font.family": "serif", "font.size": 9, "axes.linewidth": 0.7,
                     "lines.linewidth": 1.2, "lines.markersize": 4, "legend.fontsize": 8,
                     "legend.frameon": False})
SOLVERS = ["interFoam", "interIsoFoam", "interIsoGoam"]
COLOUR = {"interFoam": "#2a78d6", "interIsoFoam": "#eb6834", "interIsoGoam": "#1baf7a"}
MARKER = {"interFoam": "o", "interIsoFoam": "s", "interIsoGoam": "D"}
MESHES = [("95", "95$\\times$20"), ("190", "190$\\times$40"), ("380", "380$\\times$80"),
          ("760", "760$\\times$160")]


def load(f):
    rows = []
    for r in csv.DictReader(open(os.path.join(PD, f))):
        if r["application"] in SOLVERS and r["diverged"] != "True":
            rows.append(r)
    return rows


fixed, limit = load("runs_fixedDt.csv"), load("runs_maxCo.csv")


def series(rows, app, nx, key):
    pts = sorted((float(r["avgMaxCo"]), float(r[key])) for r in rows
                 if r["application"] == app and r["nx"] == nx and r[key] not in ("", "None"))
    return [p[0] for p in pts], [p[1] for p in pts]


fig, axes = plt.subplots(2, 4, figsize=(7.0, 4.4), sharex=True, sharey="row")
for j, (nx, label) in enumerate(MESHES):
    for i, (key, ref) in enumerate([("dHfit", 0.0), ("Umed", 1.0)]):
        ax = axes[i, j]
        for app in SOLVERS:
            x, y = series(fixed, app, nx, key)
            ax.plot(x, y, "-", color=COLOUR[app], marker=MARKER[app])
            x, y = series(limit, app, nx, key)
            ax.plot(x, y, "--", color=COLOUR[app], marker=MARKER[app], mfc="white")
        ax.axhline(ref, color="#1f1f1e", lw=0.6)
        ax.set_xscale("log")
        ax.set_xticks([0.05, 0.1, 0.2, 0.5])
        ax.set_xticklabels(["0.05", "0.1", "0.2", "0.5"])
        ax.minorticks_off()
        ax.grid(True, lw=0.4, alpha=0.5)
        if i == 0:
            ax.set_title(label, fontsize=9)
axes[0, 0].set_ylabel(r"$\Delta H$ after 15 periods [%]")
axes[1, 0].set_ylabel(r"median $\max|\mathbf{U}|/u_{crest}$")
axes[0, 0].set_ylim(-60, 15)
fig.supxlabel("time-averaged maximum Courant number", fontsize=9, y=0.01)
h = [plt.Line2D([], [], color=COLOUR[a], marker=MARKER[a], label=a) for a in SOLVERS]
h += [plt.Line2D([], [], color="#6b6b6b", ls="-", marker="o", label=r"fixed $\Delta t$"),
      plt.Line2D([], [], color="#6b6b6b", ls="--", marker="o", mfc="white", label="Courant limit")]
fig.legend(handles=h, loc="upper center", ncol=5, bbox_to_anchor=(0.5, 1.0))
fig.tight_layout(rect=(0, 0, 1, 0.94))
os.makedirs(os.path.join(HERE, "figures"), exist_ok=True)
fig.savefig(os.path.join(HERE, "figures", "waveResults.pdf"))
fig.savefig(os.path.join(HERE, "figures", "waveResults.png"), dpi=150)

# Cost under a Courant limit
MESHCOL = {"95": "#c4c4c0", "190": "#909090", "380": "#5c5c5a", "760": "#1f1f1e"}
by = {(r["application"], r["nx"], float(r["maxCo"])): r for r in limit}
fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.7), sharey=True)
for ax, key, title in [(axes[0], "steps", "(a) number of time steps"),
                       (axes[1], "clockTime", "(b) wall-clock time")]:
    for nx, label in MESHES:
        cos = sorted({c for (a, n, c) in by if n == nx}, reverse=True)
        for other, ls, mfc in [("interFoam", "-", None), ("interIsoFoam", "--", "white")]:
            y = [100*(float(by[("interIsoGoam", nx, c)][key])/float(by[(other, nx, c)][key]) - 1)
                 for c in cos]
            ax.plot(cos, y, ls, color=MESHCOL[nx], marker="o", mfc=mfc or MESHCOL[nx],
                    label=label if other == "interFoam" else None)
    ax.axhline(0, color="#1f1f1e", lw=0.6)
    ax.set_xscale("log")
    ax.set_xticks([0.05, 0.1, 0.2, 0.5])
    ax.set_xticklabels(["0.05", "0.1", "0.2", "0.5"])
    ax.minorticks_off()
    ax.grid(True, lw=0.4, alpha=0.5)
    ax.set_xlabel("maximum Courant number (limit)")
    ax.set_title(title, fontsize=9, loc="left")
axes[0].set_ylabel("change, interIsoGoam vs other [%]")
h = [plt.Line2D([], [], color=MESHCOL[n], marker="o", label=l) for n, l in MESHES]
h += [plt.Line2D([], [], color="#6b6b6b", ls="-", marker="o", label="vs interFoam"),
      plt.Line2D([], [], color="#6b6b6b", ls="--", marker="o", mfc="white", label="vs interIsoFoam")]
fig.legend(handles=h, loc="upper center", ncol=6, bbox_to_anchor=(0.5, 1.0))
fig.tight_layout(rect=(0, 0, 1, 0.9))
fig.savefig(os.path.join(HERE, "figures", "waveCost.pdf"))
fig.savefig(os.path.join(HERE, "figures", "waveCost.png"), dpi=150)
