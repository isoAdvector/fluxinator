#!/usr/bin/env python3
#
# Prints where the maximum velocity occurs for water level z_w = 0.49: after the
# first step for the basic gravity term with uncorrected schemes (interFoam), and
# after 10 steps for the geometric model with a hydrostatic start (nLayers 1 and 2).
# Gives the cell centre, the volume fraction and the ring level of the cell in the
# corrected snGrad(rho) reach field written by postProcess -time 0 -dict
# ../common/snGradRhoReach in waterLevel_10steps/g1_zero_z0.4900 (-1: no adjacent face with non-zero
# snGrad(rho)). Needs the cell centres written by postProcess -func writeCellCentres
# in waterLevel_10steps/g1_zero_z0.4900 (Allrun).
#
# Usage (from this directory):
#   python3 maxUCell.py

import re

import numpy as np

BASE = "waterLevel_10steps/g1_zero_z0.4900/0"
CASES = [("interFoam, basic, uncorrected, t = 0.01", "waterLevel_10steps/iF_unc_z0.4900/0.01"),
         ("geometric, nLayers 1, hydrostatic, t = 0.1", "waterLevel_10steps/g1_hyd_z0.4900/0.1"),
         ("geometric, nLayers 2, hydrostatic, t = 0.1", "waterLevel_10steps/g2_hyd_z0.4900/0.1")]


def internal(path):
    with open(path) as f:
        s = f.read()
    m = re.search(r"internalField\s+nonuniform\s+List<(\w+)>\s*(\d+)\s*\(", s)
    n = int(m.group(2))
    body = s[m.end():]
    if m.group(1) == "vector":
        vals = re.findall(r"\(([^()]*)\)", body)[:n]
        return np.array([[float(x) for x in v.split()] for v in vals])
    return np.array([float(x) for x in body.split(")")[0].split()[:n]])


ring = internal(f"{BASE}/snGradReachCorrected")
C = internal(f"{BASE}/C")
alpha = internal(f"{BASE}/alpha.water")
r2 = np.where(ring == 2)[0]
print("Cells in ring 2 of the corrected snGrad(rho) reach: "
      + ", ".join(f"{i} ({C[i, 0]:.2f}, {C[i, 2]:.2f})" for i in r2))
for label, path in CASES:
    mag = np.linalg.norm(internal(f"{path}/U"), axis=1)
    i = int(mag.argmax())
    print(f"{label}: max|U| = {mag[i]:.4e} m/s in cell {i} at (x, z) = ({C[i, 0]:.3f}, {C[i, 2]:.3f}),"
          f" alpha = {alpha[i]:g}, ring level {ring[i]:g}")
