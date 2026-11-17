#!/usr/bin/env python3
# Copyright (C) 2026 Johan Roenby
# SPDX-License-Identifier: GPL-3.0-or-later
"""Largest velocity of the standing wave in the air (alpha < 0.01), in the interface cells and
in the water (alpha > 0.99): after the first time step, and its peak and mean over the two
periods. Reads the fields of every time step in <dir>/<mesh>_<solver> (airWater.sh) and writes
<dir>/airWater.csv. Usage: python3 airWater.py <dir>"""
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "common"))
from foamio import read_vol  # noqa: E402

out = sys.argv[1]
rows = ["mesh,solver,first_air,first_interface,first_water,peak_air,peak_interface,peak_water,"
        "mean_air,mean_interface,mean_water"]
for mesh in ["uniform", "refined"]:
    for app in ["interFoam", "interIsoFoam", "interIsoGoam"]:
        d = os.path.join(out, f"{mesh}_{app}")
        times = sorted((t for t in os.listdir(d) if re.fullmatch(r"[0-9.e-]+", t) and float(t) > 0),
                       key=float)
        m = []
        for t in times:
            U = np.linalg.norm(read_vol(os.path.join(d, t, "U")), axis=1)
            a = read_vol(os.path.join(d, t, "alpha.water"))
            m.append([U[sel].max(initial=0) for sel in (a < 0.01, (a >= 0.01) & (a <= 0.99), a > 0.99)])
        m = np.array(m)
        rows.append(",".join([mesh, app] + [f"{v:.4g}" for v in (*m[0], *m.max(axis=0), *m.mean(axis=0))]))
open(os.path.join(out, "airWater.csv"), "w").write("\n".join(rows) + "\n")
print("\n".join(rows))
