# Copyright (C) 2026 Johan Roenby
# SPDX-License-Identifier: GPL-3.0-or-later
"""Writes 0/U of the translating Taylor-Green vortex u = U0 + sin(x - U0 t) cos z,
w = -cos(x - U0 t) sin z at t = 0 on the 64x1x64 blockMesh of [0, 2 pi]^2 (x fastest, then z).
U0 is read from the file U0 in the case. Usage: python3 writeU.py <case>"""
import os
import sys
import numpy as np

case = sys.argv[1]
U0 = float(open(os.path.join(case, "U0")).read())
n = 64
c = (np.arange(n) + 0.5)*2*np.pi/n
X, Z = np.meshgrid(c, c)                     # Z varies over rows: index = i + n*k
u = (U0 + np.sin(X)*np.cos(Z)).ravel()
w = (-np.cos(X)*np.sin(Z)).ravel()
with open(os.path.join(case, "0", "U"), "w") as f:
    f.write("FoamFile { version 2.0; format ascii; class volVectorField; object U; }\n")
    f.write("dimensions [0 1 -1 0 0 0 0];\n")
    f.write(f"internalField nonuniform List<vector> {n*n}\n(\n")
    f.writelines(f"({a:.15g} 0 {b:.15g})\n" for a, b in zip(u, w))
    f.write(");\n")
    f.write('boundaryField { ".*" { type cyclic; } frontAndBack { type empty; } }\n')
print(f"U0 = {U0}, max|U| = {np.hypot(u, w).max():.3f}")
