#!/usr/bin/env python
#
# Generates the figures for the periodic wave section of the interIsoGoam paper
# from the fixed-dt and fixed-maxCo scans. Metrics are read from the CSV files
# written by analyseScan.py (fixedDtMetrics.csv, maxCoMetrics.csv); time series are read from
# the case directories.
#
# Usage:
#   ./paperFigures.py            # writes paperFigures/*.pdf and *.png

import csv
import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from plotCrestHeight import find_extrema, wave_heights  # noqa: E402

OUT = os.path.join(HERE, "paperFigures")
FIXED_DIR = os.path.join(HERE, "gravityModelScan_fixedDt")
MAXCO_DIR = os.path.join(HERE, "gravityModelScan")
FIXED_CSV = os.path.join(HERE, "fixedDtMetrics.csv")
MAXCO_CSV = os.path.join(HERE, "maxCoMetrics.csv")
REF = "0036"

H_EXACT = 0.124982           # m, from the streamFunFoam reference
U_CREST = 0.44453480326804   # m/s
SOLVERS = ["interFoam", "interIsoFoam", "interIsoGoam"]
COLOUR = {"interFoam": "#2a78d6", "interIsoFoam": "#eb6834", "interIsoGoam": "#1baf7a"}
MARKER = {"interFoam": "o", "interIsoFoam": "s", "interIsoGoam": "D"}
NZ = [20, 40, 80, 160]
NX = {20: 95, 40: 190, 80: 380, 160: 760}

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 9,
    "axes.linewidth": 0.7,
    "lines.linewidth": 1.2,
    "lines.markersize": 4,
    "legend.fontsize": 8,
    "legend.frameon": False,
    "axes.grid": True,
    "grid.linewidth": 0.4,
    "grid.alpha": 0.5,
})
FULL, HALF = 7.0, 3.4  # figure widths [in]


# ---------------------------------------------------------------- data access
def load_csv(path, key):
    rows = []
    with open(path) as f:
        for r in csv.DictReader(f):
            out = {}
            for k, v in r.items():
                if k in ("case", "application"):
                    out[k] = v
                elif k == "diverged":
                    out[k] = (v == "True")
                elif v in ("", "None", "-"):
                    out[k] = None
                else:
                    try:
                        out[k] = float(v)
                    except ValueError:
                        out[k] = v
            out["nz"] = int(out["nz"])
            out["level"] = out[key]
            rows.append(out)
    return rows


FIX = load_csv(FIXED_CSV, "deltaT")
MC = load_csv(MAXCO_CSV, "maxCo")


def runs(rows, app, nz):
    return sorted((r for r in rows if r["application"] == app and r["nz"] == nz),
                  key=lambda r: -r["level"])


def eta(scanDir, case):
    return np.loadtxt(os.path.join(scanDir, case, "postProcessing/interfaceHeight1/0/height.dat"),
                      comments="#", unpack=True, usecols=(0, 2))


def maxU(scanDir, case):
    return np.loadtxt(os.path.join(scanDir, case, "postProcessing/fieldMinMax1/0/fieldMinMax.dat"),
                      comments="#", unpack=True, usecols=(0, 2))


def app_of(scanDir, case):
    with open(os.path.join(scanDir, case, "parmScannerParms")) as f:
        return f.read().split("application:")[1].split()[0]


def save(fig, name):
    fig.savefig(os.path.join(OUT, name + ".pdf"))
    fig.savefig(os.path.join(OUT, name + ".png"), dpi=200)
    plt.close(fig)
    print("wrote", name)


def solver_legend(ax, extra=None, **kw):
    handles = [plt.Line2D([], [], color=COLOUR[s], marker=MARKER[s], label=s) for s in SOLVERS]
    if extra:
        handles += extra
    ax.legend(handles=handles, **kw)


# ---------------------------------------------------------------- figures
def fig_eta():
    """Surface elevation at the gauge, first and last 4 s, nz 80, large and small dt."""
    settings = [(["0024", "0025", "0026"], r"$\Delta t = 5$ ms ($Co \approx 0.5$)"),
                (["0030", "0031", "0032"], r"$\Delta t = 1.25$ ms ($Co \approx 0.09$)")]
    tr, hr = eta(FIXED_DIR, REF)
    fig, axes = plt.subplots(2, 2, figsize=(FULL, 4.2), sharey=True)
    for row, (cases, label) in enumerate(settings):
        for col, (t0, t1) in enumerate([(0, 4), (26, 30)]):
            ax = axes[row, col]
            m = (tr >= t0) & (tr <= t1)
            ax.plot(tr[m], hr[m]*1e3, "k--", lw=1.0, label="exact", zorder=5)
            for c in cases:
                t, h = eta(FIXED_DIR, c)
                m = (t >= t0) & (t <= t1)
                s = app_of(FIXED_DIR, c)
                ax.plot(t[m], h[m]*1e3, color=COLOUR[s], label=s)
            ax.set_xlim(t0, t1)
            if col == 0:
                ax.set_ylabel(r"$\eta$ [mm]")
            if row == 1:
                ax.set_xlabel(r"$t$ [s]")
            ax.set_title(f"({'abcd'[2*row+col]}) {label}, " + ("first" if col == 0 else "last") + " 4 s",
                         fontsize=9, loc="left")
    handles, labels = axes[0, 1].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=4, bbox_to_anchor=(0.5, 1.0))
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    save(fig, "eta_nz80")


def fig_dH_vs_Co():
    """Trend wave height error against measured Courant number, both campaigns, per mesh."""
    fig, axes = plt.subplots(1, 4, figsize=(FULL, 2.6), sharey=True)
    for ax, nz in zip(axes, NZ):
        for s in SOLVERS:
            for rows, style, fill in ((FIX, "-", True), (MC, "--", False)):
                rr = [r for r in runs(rows, s, nz) if r["dHfit"] is not None]
                rr.sort(key=lambda r: r["avgMaxCo"])
                x = [r["avgMaxCo"] for r in rr]
                y = [r["dHfit"] for r in rr]
                e = [r["Hscatter"] for r in rr]
                ax.errorbar(x, y, yerr=e, ls=style, color=COLOUR[s], marker=MARKER[s],
                            mfc=COLOUR[s] if fill else "white", capsize=0, elinewidth=0.8)
        ax.axhline(0, color="k", lw=0.7)
        ax.set_xscale("log")
        ax.set_xlim(0.025, 0.8)
        ax.set_title(f"{NX[nz]} × {nz}", fontsize=9)
        ax.set_xlabel(r"avg. max $Co$")
    axes[0].set_ylabel(r"$\Delta H$ trend [%]")
    axes[0].set_ylim(-60, 15)
    handles = [plt.Line2D([], [], color=COLOUR[s], marker=MARKER[s], label=s) for s in SOLVERS]
    handles += [plt.Line2D([], [], color="k", ls="-", marker="o", label=r"fixed $\Delta t$"),
                plt.Line2D([], [], color="k", ls="--", marker="o", mfc="white", label=r"fixed max$Co$")]
    fig.legend(handles=handles, loc="upper center", ncol=5, bbox_to_anchor=(0.5, 1.0))
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    save(fig, "dH_vs_Co")


def fig_convergence():
    """Trend wave height error against mesh at the smallest dt / maxCo."""
    fig, ax = plt.subplots(figsize=(HALF, 2.6))
    for s in SOLVERS:
        for rows, style, fill, dx in ((FIX, "-", True, -0.03), (MC, "--", False, 0.03)):
            rr = [runs(rows, s, nz)[-1] for nz in NZ]
            x = [nz*2**dx for nz in NZ]
            ax.errorbar(x, [r["dHfit"] for r in rr], yerr=[r["Hscatter"] for r in rr], ls=style,
                        color=COLOUR[s], marker=MARKER[s], mfc=COLOUR[s] if fill else "white",
                        capsize=0, elinewidth=0.8)
    ax.axhline(0, color="k", lw=0.7)
    ax.set_xscale("log", base=2)
    ax.set_xticks(NZ)
    ax.set_xticklabels([str(n) for n in NZ])
    ax.minorticks_off()
    ax.set_xlabel(r"$n_z$ (cells over the domain height)")
    ax.set_ylabel(r"$\Delta H$ trend [%]")
    extra = [plt.Line2D([], [], color="k", ls="-", marker="o", label=r"$\Delta t_0/8$"),
             plt.Line2D([], [], color="k", ls="--", marker="o", mfc="white", label=r"max$Co$ 0.05")]
    solver_legend(ax, extra, loc="lower right", fontsize=7)
    fig.tight_layout()
    save(fig, "convergence")


def fig_crests():
    """All crest profiles of the finest runs, aligned on crest time."""
    cases = ["0046", "0047", "0048"]
    tr, hr = eta(FIXED_DIR, REF)
    rc = find_extrema(tr, hr, 1, 0.1*(hr.max()-hr.min()))
    tw = np.linspace(-0.3, 0.3, 241)
    ref_profile = np.interp(rc[7, 0] + tw, tr, hr)*1e3
    fig, axes = plt.subplots(1, 3, figsize=(FULL, 2.4), sharey=True)
    for ax, c in zip(axes, cases):
        s = app_of(FIXED_DIR, c)
        t, h = eta(FIXED_DIR, c)
        cr = find_extrema(t, h, 1, 0.1*(h.max()-h.min()))
        n = len(cr)
        for j, (tc, _) in enumerate(cr[1:], start=1):
            if tc + 0.3 > t[-1]:
                continue
            ax.plot(tw, np.interp(tc + tw, t, h)*1e3, color=COLOUR[s], lw=0.9, alpha=0.3 + 0.7*j/n)
        ax.plot(tw, ref_profile, "k--", lw=1.0)
        r = next(x for x in FIX if x["case"] == c)
        ax.set_title(f"{s} ({r['crestPeaks']:.1f} peaks/crest)", fontsize=9)
        ax.set_xlabel(r"$t - t_{crest}$ [s]")
        ax.set_xlim(-0.3, 0.3)
    axes[0].set_ylabel(r"$\eta$ [mm]")
    axes[0].set_ylim(55, 92)
    fig.tight_layout()
    save(fig, "crests_nz160")


def fig_perwave():
    """Per-wave height error with quadratic trend (nz 160, dt0/4), and crest timing (nz 80, dt0)."""
    tr, hr = eta(FIXED_DIR, REF)
    rcr = find_extrema(tr, hr, 1, 0.1*(hr.max()-hr.min()))[:, 0]
    fig, axes = plt.subplots(1, 2, figsize=(FULL, 2.5))
    for c in ["0043", "0044", "0045"]:
        s = app_of(FIXED_DIR, c)
        t, h = eta(FIXED_DIR, c)
        thr = 0.1*(h.max()-h.min())
        H = wave_heights(find_extrema(t, h, 1, thr), find_extrema(t, h, -1, thr))
        rel = (H[:, 1] - H_EXACT)/H_EXACT*100
        p = np.polyfit(H[:, 0], rel, 2)
        tt = np.linspace(H[0, 0], H[-1, 0], 100)
        axes[0].plot(H[:, 0], rel, marker=MARKER[s], color=COLOUR[s], lw=0.7, label=s)
        axes[0].plot(tt, np.polyval(p, tt), "--", color=COLOUR[s], lw=1.0)
    axes[0].axhline(0, color="k", lw=0.7)
    axes[0].set_xlabel(r"$t$ [s]")
    axes[0].set_ylabel(r"$(H_n - H)/H$ [%]")
    axes[0].set_title(r"(a) $760\times160$, $\Delta t_0/4$: per wave and quadratic trend", fontsize=9, loc="left")
    axes[0].legend(loc="lower left")
    for c in ["0024", "0025", "0026"]:
        s = app_of(FIXED_DIR, c)
        t, h = eta(FIXED_DIR, c)
        cr = find_extrema(t, h, 1, 0.1*(h.max()-h.min()))[:, 0]
        # Match each crest to the nearest exact crest, but only within half a
        # period; the exact record may end before the last simulated crest
        pairs = [(tc, tc - rcr[np.argmin(abs(rcr - tc))]) for tc in cr]
        pairs = [(tc, e) for tc, e in pairs if abs(e) < 1.0]
        axes[1].plot([p[0] for p in pairs], [p[1] for p in pairs], marker=MARKER[s],
                     color=COLOUR[s], lw=0.7, label=s)
    axes[1].axhline(0, color="k", lw=0.7)
    axes[1].set_xlabel(r"$t$ [s]")
    axes[1].set_ylabel(r"$t_{crest} - t_{crest,exact}$ [s]")
    axes[1].set_title(r"(b) $380\times80$, $\Delta t_0$: crest arrival time", fontsize=9, loc="left")
    fig.tight_layout()
    save(fig, "perwave_phase")


def binned_max(t, v, width=0.1, tEnd=30.0):
    edges = np.arange(0, tEnd + 1e-9, width)
    idx = np.digitize(t, edges) - 1
    out = np.full(len(edges) - 1, np.nan)
    for i in range(len(edges) - 1):
        sel = v[idx == i]
        if len(sel):
            out[i] = sel.max()
    return edges[:-1] + width/2, out


def fig_umax():
    """(a) max|U| over time, nz 160, dt0/8. (b) median max|U| against Courant number, all runs."""
    fig, axes = plt.subplots(1, 2, figsize=(FULL, 2.7))
    tr, ur = maxU(FIXED_DIR, REF)
    x, y = binned_max(tr, ur/U_CREST)
    axes[0].plot(x, y, "k--", lw=1.0, label="exact")
    for c in ["0046", "0047", "0048"]:
        s = app_of(FIXED_DIR, c)
        t, u = maxU(FIXED_DIR, c)
        x, y = binned_max(t, u/U_CREST)
        axes[0].plot(x, y, color=COLOUR[s], lw=0.9, label=s)
    axes[0].set_xlim(0, 30)
    axes[0].set_ylim(0.9, 3.6)
    axes[0].set_xlabel(r"$t$ [s]")
    axes[0].set_ylabel(r"max$|\mathbf{U}|/u_{crest}$ (max per 0.1 s)")
    axes[0].set_title(r"(a) $760\times160$, $\Delta t_0/8$", fontsize=9, loc="left")
    axes[0].legend(loc="upper left", ncol=2)
    for s in SOLVERS:
        for rows, fill in ((FIX, True), (MC, False)):
            rr = [r for r in rows if r["application"] == s and r["Umed"] is not None]
            axes[1].plot([r["avgMaxCo"] for r in rr], [r["Umed"] for r in rr], ls="none",
                         marker=MARKER[s], color=COLOUR[s], mfc=COLOUR[s] if fill else "white")
    axes[1].axhline(1, color="k", lw=0.7)
    axes[1].set_xscale("log")
    axes[1].set_xlabel(r"avg. max $Co$")
    axes[1].set_ylabel(r"median max$|\mathbf{U}|/u_{crest}$")
    axes[1].set_title("(b) all runs, all meshes", fontsize=9, loc="left")
    extra = [plt.Line2D([], [], color="k", ls="none", marker="o", label=r"fixed $\Delta t$"),
             plt.Line2D([], [], color="k", ls="none", marker="o", mfc="white", label=r"fixed max$Co$")]
    solver_legend(axes[1], extra, loc="upper left", fontsize=7)
    fig.tight_layout()
    save(fig, "spurious_velocity")


def fig_dt():
    """Adaptive time step over time, nz 160, maxCo 0.05."""
    fig, ax = plt.subplots(figsize=(HALF, 2.4))
    edges = np.arange(0, 30 + 1e-9, 0.1)
    for c in ["0046", "0047", "0048"]:
        s = app_of(MAXCO_DIR, c)
        t, _ = eta(MAXCO_DIR, c)
        counts = np.bincount(np.digitize(t, edges) - 1, minlength=len(edges) - 1)[:len(edges) - 1]
        dt = np.where(counts > 0, 100.0/np.maximum(counts, 1), np.nan)
        ax.plot(edges[:-1][3:] + 0.05, dt[3:], color=COLOUR[s], lw=0.9, label=s)
    ax.set_xlim(0, 30)
    ax.set_xlabel(r"$t$ [s]")
    ax.set_ylabel(r"mean $\Delta t$ per 0.1 s [ms]")
    ax.legend(loc="lower left")
    fig.tight_layout()
    save(fig, "dt_maxCo")


def interp_fixed(r):
    pts = sorted((math.log(q["deltaT"]), q["dHfit"]) for q in FIX
                 if q["application"] == r["application"] and q["nz"] == r["nz"] and q["dHfit"] is not None)
    x = math.log(30.0/r["steps"])
    if x < pts[0][0] or x > pts[-1][0]:
        return None
    for (x0, y0), (x1, y1) in zip(pts, pts[1:]):
        if x0 <= x <= x1:
            return y0 + (y1 - y0)*(x - x0)/(x1 - x0)
    return None


def fig_agreement():
    fig, ax = plt.subplots(figsize=(HALF, 3.2))
    diffs = []
    for s in SOLVERS:
        for r in (q for q in MC if q["application"] == s):
            yi = interp_fixed(r)
            if yi is None:
                continue
            diffs.append(abs(r["dHfit"] - yi))
            ax.errorbar(yi, r["dHfit"], yerr=r["Hscatter"], marker=MARKER[s], color=COLOUR[s],
                        ls="none", elinewidth=0.8)
    ax.plot([-60, 15], [-60, 15], "k-", lw=0.7)
    ax.set_xlim(-60, 15)
    ax.set_ylim(-60, 15)
    ax.set_aspect("equal")
    ax.set_xlabel(r"$\Delta H$ trend, fixed $\Delta t$ interpolated [%]")
    ax.set_ylabel(r"$\Delta H$ trend, fixed max$Co$ [%]")
    solver_legend(ax, loc="upper left")
    fig.tight_layout()
    save(fig, "agreement")
    d = sorted(diffs)
    print(f"agreement: n={len(d)} median={d[len(d)//2]:.2f} within2={sum(x <= 2 for x in d)}")


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    fig_eta()
    fig_dH_vs_Co()
    fig_convergence()
    fig_crests()
    fig_perwave()
    fig_umax()
    fig_dt()
    fig_agreement()
