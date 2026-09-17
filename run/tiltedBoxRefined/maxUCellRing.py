#!/usr/bin/env python3
#
# Prints where the maximum velocity occurs after 10 time steps for the geometric
# gravity model from a hydrostatic start (nLayers 1 and 2), and the ring level
# of that cell in the corrected snGrad(rho) reach field written by
# postProcess -time 0 -dict ../snGradRhoReach in zeroStartScan/0000
# (-1: no adjacent face with non-zero snGrad(rho)).
#
# Usage (from this directory):
#   python3 maxUCellRing.py

import re

import numpy as np

REACH = "zeroStartScan/0000/0/snGradReachCorrected"
CASES = [("nLayers 1", "nLayersHydrostaticScan/0000/0.1/U"),
         ("nLayers 2", "nLayersHydrostaticScan/0001/0.1/U")]


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


ring = internal(REACH)
print(f"Cells in ring 2 of the corrected snGrad(rho) reach: {np.where(ring == 2)[0].tolist()}")
for label, path in CASES:
    mag = np.linalg.norm(internal(path), axis=1)
    i = int(mag.argmax())
    print(f"{label}: max|U| = {mag[i]:.4e} m/s in cell {i}, ring level {ring[i]:g}")
