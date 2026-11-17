#!/usr/bin/env python3
# Copyright (C) 2026 Johan Roenby
# SPDX-License-Identifier: GPL-3.0-or-later
"""Part I, Table tab:standingWave: standing wave at the still water level z = 0.5 m on the
uniform and the refined mesh: max|U| after one time step, its peak over the first twenty
steps, its peak and mean over the first two wave periods, and the number of cells with
1e-6 < alpha < 1 - 1e-6 at t = 1.2 s, about one period.
Inputs (run/standingWave/paperScan.sh): <dir>/history[Uniform]_<solver>.csv and
<dir>/hist[Uniform]_<solver>/1.2/alpha.water, with <dir> = solvers_plicRDF (the table) and
solvers (isoAlpha: largest relative change of each velocity, printed below the table).
Output: tab_standingWave.txt. Usage: python3 make.py"""
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
from common import case, write_table  # noqa: E402

SOLVERS = ["interFoam", "interIsoFoam", "interIsoGoam"]
TAG = {"uniform": "Uniform", "refined": ""}


def alpha_band(path):
    s = open(path).read()
    m = re.search(r"internalField\s+nonuniform\s+List<scalar>\s*(\d+)\s*\(", s)
    a = np.array(s[m.end():].split()[:int(m.group(1))], float)
    return int(((a > 1e-6) & (a < 1 - 1e-6)).sum())


def velocities(d, mesh, app):
    """max|U| after step 1, peak over 20 steps, peak and mean over the two periods"""
    t, u = np.loadtxt(case("standingWave", d, f"history{TAG[mesh]}_{app}.csv"),
                      delimiter=",", unpack=True)
    return np.array([u[0], u[:20].max(), u.max(), u.mean()])


out = ["Standing wave at z = 0.5 m (plicRDF): max|U| [m/s] and interface band",
       f"{'mesh':8s} {'solver':13s} {'1 step':>9s} {'peak20':>9s} {'2T peak':>9s} {'2T mean':>9s} {'band':>6s}"]
change = np.zeros(4)
for mesh in TAG:
    for app in SOLVERS:
        v = velocities("solvers_plicRDF", mesh, app)
        band = alpha_band(case("standingWave", "solvers_plicRDF", f"hist{TAG[mesh]}_{app}", "1.2", "alpha.water"))
        out.append(f"{mesh:8s} {app:13s} {v[0]:9.3f} {v[1]:9.3f} {v[2]:9.2f} {v[3]:9.2f} {band:6d}")
        if app != "interFoam":
            change = np.maximum(change, np.abs(velocities("solvers", mesh, app)/v - 1))
out += ["", "isoAlpha instead of plicRDF (interIsoFoam, interIsoGoam): largest relative change [%]",
        "  1 step %.0f, peak20 %.0f, 2T peak %.0f, 2T mean %.0f" % tuple(100*change)]
out += ["", "max|U| [m/s] in the air (alpha < 0.01), the interface cells and the water (alpha > 0.99)",
        "(run/standingWave/airWater.sh: the same runs with every time step written)",
        f"{'mesh':8s} {'solver':13s} {'1 step: air':>12s} {'interface':>10s} {'water':>8s}"
        f" {'2T peak: air':>13s} {'interface':>10s} {'water':>8s}"]
for r in np.genfromtxt(case("standingWave", "solvers_airWater", "airWater.csv"), delimiter=",",
                       names=True, dtype=None, encoding=None):
    out.append(f"{r['mesh']:8s} {r['solver']:13s} {r['first_air']:12.3f} {r['first_interface']:10.3f}"
               f" {r['first_water']:8.4f} {r['peak_air']:13.2f} {r['peak_interface']:10.2f} {r['peak_water']:8.3f}")
write_table(HERE, "tab_standingWave", out)
