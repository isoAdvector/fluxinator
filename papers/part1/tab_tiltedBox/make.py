#!/usr/bin/env python3
# Copyright (C) 2026 Johan Roenby
# SPDX-License-Identifier: GPL-3.0-or-later
"""Part I, Table tab:tiltedBox: maximum velocity magnitude [m/s] in the tilted box
at t = 0.01 s, at t = 10 s and its peak over t <= 10 s.
Inputs: run/tiltedBox/gravityModels_1000steps (cases 0000-0002) and
run/tiltedBox/nLayers_1000steps (geometric with nLayers 0, 1, 2).
Output: tab_tiltedBox.txt. Usage: python3 make.py"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
from common import case  # noqa: E402

ROWS = [("interFoam, basic", "gravityModels_1000steps", "0000"),
        ("interIsoFoam, basic", "gravityModels_1000steps", "0001"),
        ("interIsoGoam, constantHeight", "gravityModels_1000steps", "0002"),
        ("interIsoGoam, geometric, nL = 0", "nLayers_1000steps", "0000"),
        ("interIsoGoam, geometric, nL = 1", "nLayers_1000steps", "0001"),
        ("interIsoGoam, geometric, nL = 2", "nLayers_1000steps", "0002")]

out = ["Maximum velocity magnitude [m/s] in the tilted box",
       f"{'configuration':34s} {'t = 0.01 s':>11s} {'t = 10 s':>11s} {'peak':>11s}   source"]
for label, study, c in ROWS:
    f = case("tiltedBox", study, c, "postProcessing", "fieldMinMax1", "0", "fieldMinMax.dat")
    t, u = np.loadtxt(f, comments="#", usecols=(0, 2), unpack=True, ndmin=2)
    vals = [u[np.argmin(abs(t - T))] for T in (0.01, 10)] + [u.max()]
    out.append(f"{label:34s}" + "".join(f" {v:11.3e}" for v in vals) + f"   run/tiltedBox/{study}/{c}")
open(os.path.join(HERE, "tab_tiltedBox.txt"), "w").write("\n".join(out) + "\n")
print("\n".join(out))
