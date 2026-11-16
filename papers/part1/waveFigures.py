#!/usr/bin/env python
# Copyright (C) 2026 Johan Roenby
# SPDX-License-Identifier: GPL-3.0-or-later
#
# Shared code of the Part I progressive-wave figures fig_eta, fig_crests and
# fig_spurious (each make.py sets OUT and calls one fig_* function). Metrics are
# read from run/progressiveWave/data/solvers_{fixedDt,maxCo}.csv (analyseScan.py);
# time series from the case directories of run/progressiveWave/solvers_*.

import csv
import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
WAVE = os.path.join(HERE, "..", "..", "run", "progressiveWave")
sys.path.insert(0, WAVE)
from plotCrestHeight import find_extrema, wave_heights  # noqa: E402
from common import paper_style  # noqa: E402

OUT = HERE  # set by the calling make.py
FIXED_DIR = os.path.join(WAVE, "solvers_fixedDt")
MAXCO_DIR = os.path.join(WAVE, "solvers_maxCo")
FIXED_CSV = os.path.join(WAVE, "data", "solvers_fixedDt.csv")
MAXCO_CSV = os.path.join(WAVE, "data", "solvers_maxCo.csv")
REF = "0036"

H_EXACT = 0.124982           # m, from the streamFunFoam reference
U_CREST = 0.44453480326804   # m/s
WAVELENGTH = 3.7873860376294 # m, the domain length (one wavelength, cyclic)
PERIOD = 2.0                 # s
C_EXACT = WAVELENGTH/PERIOD  # m/s, the exact celerity
SOLVERS = ["interFoam", "interIsoFoam", "interIsoGoam"]
COLOUR = {"interFoam": "#DDAA33", "interIsoFoam": "#BB5566", "interIsoGoam": "#004488"}  # Paul Tol high-contrast
MARKER = {"interFoam": "o", "interIsoFoam": "s", "interIsoGoam": "D"}
NZ = [20, 40, 80, 160]
NX = {20: 95, 40: 190, 80: 380, 160: 760}

plt.rcParams.update({
    "axes.linewidth": 0.7,
    "lines.linewidth": 1.2,
    "lines.markersize": 4,
    "legend.frameon": False,
    "axes.grid": True,
    "grid.linewidth": 0.4,
    "grid.alpha": 0.5,
})
FULL, HALF = 7.0, 3.4  # figure widths [in]
paper_style(scale=5.0/FULL)  # FULL is printed at the text width


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
    return np.loadtxt(os.path.join(scanDir, case, "postProcessing/surfacePosition1/0/height.dat"),
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
    # Zoom on the crest near t = 28 s in the last 4 s ([t0, t1, eta0, eta1], s and mm),
    # drawn above the trough at t = 27 s, where all curves stay below 10 mm
    zoom = [[27.82, 28.06, 66, 75], [27.95, 28.28, 72, 86]]
    ylim = (-50, 106)  # mm, fixed so that the inset can be placed in data units
    insetBounds = [26.55, 22, 1.05, 76]  # t, eta, width, height (s, mm)
    insetBounds = [(insetBounds[0] - 26)/4, (insetBounds[1] - ylim[0])/(ylim[1] - ylim[0]),
                   insetBounds[2]/4, insetBounds[3]/(ylim[1] - ylim[0])]  # axes fraction
    tr, hr = eta(FIXED_DIR, REF)
    fig, axes = plt.subplots(2, 2, figsize=(FULL, 4.2), sharey=True)
    for row, (cases, label) in enumerate(settings):
        for col, (t0, t1) in enumerate([(0, 4), (26, 30)]):
            ax = axes[row, col]
            axs = [ax]
            if col == 1:
                axs.append(ax.inset_axes(insetBounds))
            m = (tr >= t0) & (tr <= t1)
            for a in axs:
                a.plot(tr[m], hr[m]*1e3, "k--", lw=1.0, label="exact", zorder=5)
            for c in cases:
                t, h = eta(FIXED_DIR, c)
                m = (t >= t0) & (t <= t1)
                s = app_of(FIXED_DIR, c)
                for a in axs:
                    a.plot(t[m], h[m]*1e3, color=COLOUR[s], label=s)
            if col == 1:
                inset = axs[1]
                inset.set_xlim(*zoom[row][:2])
                inset.set_ylim(*zoom[row][2:])
                inset.set_xticks([])
                inset.set_yticks([])
                inset.grid(False)
                ax.indicate_inset_zoom(inset, edgecolor="0.4")
            ax.set_xlim(t0, t1)
            ax.set_ylim(*ylim)
            if col == 0:
                ax.set_ylabel(r"$\eta$ [mm]")
            if row == 1:
                ax.set_xlabel(r"$t$ [s]")
            ax.set_title(f"({'abcd'[2*row+col]}) {label}, " + ("first" if col == 0 else "last") + " 4 s",
                         loc="left")
    handles, labels = axes[0, 1].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=4, bbox_to_anchor=(0.5, 1.0))
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    save(fig, "eta")


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
        ax.set_title(f"{s}: {r['crestPeaks']:.1f} peaks")
        ax.set_xlabel(r"$t - t_{crest}$ [s]")
        ax.set_xlim(-0.3, 0.3)
    axes[0].set_ylabel(r"$\eta$ [mm]")
    axes[0].set_ylim(55, 92)
    fig.tight_layout()
    save(fig, "crests")


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
    """max|U| over time, nz 80, dt0/8 (a mesh and time step with the typical ordering of the
    solvers; on the 760x160 mesh at the smallest time steps interFoam is higher). The median
    max|U| of all runs is shown in fig_waveResults."""
    fig, ax = plt.subplots(figsize=(0.65*FULL, 2.7))
    tr, ur = maxU(FIXED_DIR, REF)
    x, y = binned_max(tr, ur/U_CREST)
    ax.plot(x, y, "k--", lw=1.0, label="exact")
    for c in ["0033", "0034", "0035"]:
        s = app_of(FIXED_DIR, c)
        t, u = maxU(FIXED_DIR, c)
        x, y = binned_max(t, u/U_CREST)
        ax.plot(x, y, color=COLOUR[s], lw=0.9, label=s)
    ax.set_xlim(0, 30)
    ax.set_ylim(0.9, 2.2)
    ax.set_xlabel(r"$t$ [s]")
    ax.set_ylabel(r"max$|\mathbf{U}|/u_{crest}$ (max per 0.1 s)")
    ax.set_title(r"$380\times80$, $\Delta t_0/8$", loc="left")
    ax.legend(loc="upper left", ncol=2)
    fig.tight_layout()
    save(fig, "spurious")


