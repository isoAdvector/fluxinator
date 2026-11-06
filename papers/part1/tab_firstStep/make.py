#!/usr/bin/env python3
# Copyright (C) 2026 Johan Roenby
# SPDX-License-Identifier: GPL-3.0-or-later
"""Part I, Table tab:firstStep: first-step spurious acceleration max|U|/(dt g) of the
standard (basic) gravity term in the tilted box, per mesh, min ... max over
dt = 1e-4 ... 1 s. Input: run/tiltedBox/data/meshTimeStep_1step.txt (written by
run/tiltedBox/writeData.py from run/tiltedBox/meshTimeStep_1step).
Output: tab_firstStep.txt. Usage: python3 make.py"""
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
from common import case  # noqa: E402

txt = open(case("tiltedBox", "data", "meshTimeStep_1step.txt")).read().split("\n")
rows = {}
for i, l in enumerate(txt):
    m = re.match(r"\d{4}\s+nx: (\d+).*application: (\S+) deltaT: (\S+)", l)
    if m:
        u = float(re.search(r"max\|U\| t=\S+: (\S+)", txt[i + 1]).group(1))
        rows.setdefault(int(m.group(1)), []).append((m.group(2), float(m.group(3)), u))
out = ["First time step in the tilted box, standard (basic) gravity term: max|U|/(dt g), g = 9.81 m/s^2",
       "Source: run/tiltedBox/meshTimeStep_1step (run/tiltedBox/data/meshTimeStep_1step.txt)",
       f'{"mesh":>10s} {"interFoam":>22s} {"interIsoFoam (basic)":>22s}   (min ... max over dt = 1e-4 ... 1 s)']
for n in sorted(rows):
    s = []
    for app in ("interFoam", "interIsoGoam"):
        a = [u/dt/9.81 for ap, dt, u in rows[n] if ap == app]
        s.append(f"{min(a):.4f} ... {max(a):.4f}")
    out.append(f"{n}x{n:<7d} {s[0]:>22s} {s[1]:>22s}")
open(os.path.join(HERE, "tab_firstStep.txt"), "w").write("\n".join(out) + "\n")
print("\n".join(out))
