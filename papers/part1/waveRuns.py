# Copyright (C) 2026 Johan Roenby
# SPDX-License-Identifier: GPL-3.0-or-later
"""Progressive-wave results of Part I as figures, from the metrics of all runs
(run/progressiveWave/data/solvers_fixedDt.csv and solvers_maxCo.csv, written by analyseScan.py).
fig_waveResults: for each mesh, (top) wave height error after 15 periods (dHfit) and
(bottom) median maximum velocity over t > 2 s relative to the crest velocity (Umed), against the
time-averaged maximum Courant number; filled: fixed time step, open: Courant limit.
fig_waveCost: runs under a Courant limit, change of interIsoGoam relative to interFoam
(solid) and interIsoFoam (dashed) in (a) the number of time steps and (b) the wall-clock time,
against the Courant limit, one line per mesh.
Shared data and style of papers/part1/fig_waveResults and fig_waveCost."""
import csv
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
PD = os.path.join(HERE, "..", "..", "run", "progressiveWave", "data")
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


fixed, limit = load("solvers_fixedDt.csv"), load("solvers_maxCo.csv")


def series(rows, app, nx, key):
    pts = sorted((float(r["avgMaxCo"]), float(r[key])) for r in rows
                 if r["application"] == app and r["nx"] == nx and r[key] not in ("", "None"))
    return [p[0] for p in pts], [p[1] for p in pts]


