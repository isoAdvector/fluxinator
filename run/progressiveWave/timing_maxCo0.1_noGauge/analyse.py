# Copyright (C) 2026 Johan Roenby
# SPDX-License-Identifier: GPL-3.0-or-later
"""Time per step of each timing run, over 2 <= t <= 6 s, from ExecutionTime in the log.
Writes ../data/timing_maxCo0.1.csv. Usage: python3 analyse.py"""
import csv
import glob
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
rows = []
for d in sorted(glob.glob(os.path.join(HERE, "*_r[123]"))):
    solver, mesh, rep = os.path.basename(d).split("_")
    log = [l for l in glob.glob(d + "/log.*") if not re.search("blockMesh|streamFun|dictator|Allrun", l)][0]
    t, ex = None, []
    for line in open(log):
        if line.startswith("Time = "):
            t = float(line.split()[2])
        elif line.startswith("ExecutionTime") and t is not None:
            ex.append((t, float(line.split()[2])))
    i = next(k for k, r in enumerate(ex) if r[0] >= 2.0)
    assert ex[-1][0] > 5.999, d
    nx, nz = mesh.split("x")
    rows.append([solver, nx, nz, rep[1:], f"{1e3*(ex[-1][1] - ex[i][1])/(len(ex) - 1 - i):.3f}"])
with open(os.path.join(HERE, "../data/timing_maxCo0.1.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["solver", "nx", "nz", "repeat", "msPerStep"])
    w.writerows(rows)
