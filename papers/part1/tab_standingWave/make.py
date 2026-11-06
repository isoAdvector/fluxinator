#!/usr/bin/env python3
# Copyright (C) 2026 Johan Roenby
# SPDX-License-Identifier: GPL-3.0-or-later
"""Part I, Table tab:standingWave: standing wave, max|U| after one time step and its peak
over the first twenty steps (median over twelve interface offsets), peak and mean of
max|U| over the first two wave periods on the refined mesh, and the number of cells with
1e-6 < alpha < 1 - 1e-6 at t = 1.2 s, about one period (history run, middle offset).
Inputs: run/standingWave/solvers/scan.csv, solvers/history_<solver>.csv and
solvers/hist_<solver>/1.2/alpha.water (paperScan.sh).
Output: tab_standingWave.txt. Usage: python3 make.py"""
import csv
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
from common import case, write_table  # noqa: E402

SOLVERS = ["interFoam", "interIsoFoam", "interIsoGoam"]
rows = list(csv.DictReader(open(case("standingWave", "solvers", "scan.csv"))))


def alpha_band(path):
    s = open(path).read()
    m = re.search(r"internalField\s+nonuniform\s+List<scalar>\s*(\d+)\s*\(", s)
    a = np.array(s[m.end():].split()[:int(m.group(1))], float)
    return int(((a > 1e-6) & (a < 1 - 1e-6)).sum())


out = ["Standing wave (median over interface offsets; history: refined mesh, two periods)",
       f"{'mesh':8s} {'solver':13s} {'1 step':>9s} {'peak20':>9s} {'2T peak':>9s} {'2T mean':>9s} {'band':>6s}"]
for mesh in ["uniform", "refined"]:
    for app in SOLVERS:
        r = [x for x in rows if x["mesh"] == mesh and x["app"] == app]
        first = np.median([float(x["first"]) for x in r])
        peak = np.median([float(x["peak20"]) for x in r])
        line = f"{mesh:8s} {app:13s} {first:9.3f} {peak:9.3f}"
        if mesh == "refined":
            t, u = np.loadtxt(case("standingWave", "solvers", f"history_{app}.csv"), delimiter=",", unpack=True)
            band = alpha_band(case("standingWave", "solvers", f"hist_{app}", "1.2", "alpha.water"))
            line += f" {u.max():9.2f} {u.mean():9.2f} {band:6d}"
        out.append(line)
write_table(HERE, "tab_standingWave", out)
