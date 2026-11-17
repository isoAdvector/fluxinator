# Copyright (C) 2026 Johan Roenby
# SPDX-License-Identifier: GPL-3.0-or-later
"""Kinematic test rig (applications/test/convectionTest, run/prescribedFlows): per variant and case
the maximum over all time steps of the largest velocity error [m/s] and the mean over the
time steps of the rms error in interface cells. Writes data/rig.txt (tables tab:rigUniform,
tab:rigFaces) and tab_rigConvergence/rigConvergence.txt (run by the make.py of these tables)."""
import os
import numpy as np
from common import RIG, write

V = [("implEuler_new", "implicit, Euler, rho^{n+1}"), ("implEuler_old", "implicit, Euler, rho^n"),
     ("explEuler_new", "explicit, Euler, rho^{n+1}"), ("explEuler_old", "explicit, Euler, rho^n"),
     ("implBackward_new", "implicit, backward, rho^{n+1}"),
     ("explLinear_new", "explicit, linear"), ("explLW_new", "explicit, LW (limitedLinearV)"),
     ("explIU_new", "explicit, interfaceUpwind"), ("explLWIU_new", "explicit, LW + interfaceUpwind"),
     ("implIU_new", "implicit, interfaceUpwind")]


def stats(case, v):
    f = os.path.join(RIG, case, "postProcessing", "convectionTest", v + ".dat")
    if not os.path.exists(f):
        return None
    d = np.loadtxt(f, ndmin=2)
    return d[:, 1].max(), d[:, 4].mean()


def table(cases, title):
    lines = [title, f"{'variant':34s}" + "".join(f"{c:>24s}" for c in cases)]
    for v, lab in V:
        row = []
        for c in cases:
            s = stats(c, v)
            row.append("--" if s is None else f"{s[0]:.2e} ({s[1]:.1e})")
        lines.append(f"{lab:34s}" + "".join(f"{x:>24s}" for x in row))
    return lines


cases = ["uniform", "acceleration", "accelerationCarried", "shear", "shearCarried", "rotation", "rotationCarried"]
lines = ["Kinematic rig, 50x50, dt 0.005 s, 1 s, water/air 1000/1, disc r 0.15 m.",
         "Entry: max over steps of max |U - u_exact| [m/s] (mean over steps of interface rms).",
         "Flows: uniform U0 (0.2,0,0.05); acceleration U0 (0.1,0,0.05) + a (0.5,0,0); shear u = (z-0.5, 0, -0.2), disc at (0.5, 0.75);",
         "rotation solid body Omega 2 about (0.5,0.5). Reference: u at the cell centre. Mode: one step (oneStep) or",
         "carried (carried, cases ...Carried). Sources: run/prescribedFlows/<case>/postProcessing/convectionTest/<variant>.dat", ""]
lines += table(cases[:5], "-- time scheme, time level, shear")
lines += [""] + table(cases[5:], "-- rotation")
write("rig.txt", lines)

conv = ["Mesh convergence (dt scaled with the cell size, 25: 0.01, 50: 0.005, 100: 0.0025)",
        "Entry: max over steps of max error [m/s] (mean interface rms)", ""]
for base in ["shear", "shearCarried", "rotation", "rotationCarried"]:
    conv += table([f"{base}_n25", base, f"{base}_n100"], f"-- {base}: 25x25, 50x50, 100x100") + [""]
write("rigConvergence.txt", conv)
