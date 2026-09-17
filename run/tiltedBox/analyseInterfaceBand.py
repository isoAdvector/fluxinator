#!/usr/bin/env python3
#
# Interface band and numerical vapour in the tilted box with the standard
# gravity term (interFoam and interIsoFoam), from the volume fraction fields of
# gravityModelScan_100deltaT (t = 1 s) and gravityModelScan_1000deltaT
# (t = 10 s). For each run it prints
#   - the number of cells with 1e-3 < alpha < 0.999 (interface band),
#   - for the cells whose centre lies more than 0.06 m (1.5 cell sizes) above
#     the initial interface plane z = 0.5 m: their number, the number with
#     alpha > 1e-9, their mean and maximum alpha, the resulting mean air
#     density and the density ratio across the interface.
# Cell centres are computed from constant/polyMesh as the mean of the cell
# vertices, which is exact for the parallelogram cells of this mesh.
#
# Usage (from this directory):
#   python3 analyseInterfaceBand.py

import os
import re

import numpy as np

RHO_WATER = 1000.0
RHO_AIR = 1.0
Z_INTERFACE = 0.5
HEIGHT_ABOVE = 0.06  # m, 1.5 cell sizes of the 25 x 25 mesh of the 1 m box

RUNS = [  # (label, scan, case, time)
    ("interFoam", "gravityModelScan_100deltaT", "0000", "1"),
    ("interIsoFoam", "gravityModelScan_100deltaT", "0001", "1"),
    ("interFoam", "gravityModelScan_1000deltaT", "0000", "10"),
    ("interIsoFoam", "gravityModelScan_1000deltaT", "0001", "10"),
]


def foam_list(path):
    """Body of the first ascii list in an OpenFOAM file."""
    with open(path) as f:
        s = f.read()
    s = s[s.index("}") + 1:]  # skip the FoamFile header
    m = re.search(r"\n(\d+)\s*\n\(", s)
    return int(m.group(1)), s[m.end():]


def read_points(path):
    n, body = foam_list(path)
    vals = re.findall(r"\(([^()]*)\)", body)[:n]
    return np.array([[float(x) for x in v.split()] for v in vals])


def read_faces(path):
    n, body = foam_list(path)
    return [[int(x) for x in f.split()] for f in re.findall(r"\d+\(([^()]*)\)", body)[:n]]


def read_labels(path):
    n, body = foam_list(path)
    return [int(x) for x in body.split(")")[0].split()[:n]]


def cell_centres(meshDir):
    points = read_points(os.path.join(meshDir, "points"))
    faces = read_faces(os.path.join(meshDir, "faces"))
    owner = read_labels(os.path.join(meshDir, "owner"))
    neighbour = read_labels(os.path.join(meshDir, "neighbour"))
    nCells = max(owner) + 1
    cellPoints = [set() for _ in range(nCells)]
    for i, face in enumerate(faces):
        cellPoints[owner[i]].update(face)
        if i < len(neighbour):
            cellPoints[neighbour[i]].update(face)
    return np.array([points[sorted(p)].mean(axis=0) for p in cellPoints])


def internal_scalar(path):
    with open(path) as f:
        s = f.read()
    m = re.search(r"internalField\s+nonuniform\s+List<scalar>\s*(\d+)\s*\(", s)
    if m:
        n = int(m.group(1))
        return np.array([float(x) for x in s[m.end():].split(")")[0].split()[:n]])
    m = re.search(r"internalField\s+uniform\s+(\S+);", s)
    return float(m.group(1))


def band(alpha):
    return int(((alpha > 1e-3) & (alpha < 0.999)).sum())


def main():
    print(f"Cells more than {HEIGHT_ABOVE} m above z = {Z_INTERFACE} m count as air region.")
    print(f"{'solver':13s} {'t [s]':>5s} {'band':>5s} {'nAir':>5s} {'wet':>5s} "
          f"{'mean alpha':>11s} {'max alpha':>10s} {'rho_air':>8s} {'ratio':>7s} {'reduction %':>11s}")
    for label, scan, case, time in RUNS:
        caseDir = os.path.join(scan, case)
        C = cell_centres(os.path.join(caseDir, "constant", "polyMesh"))
        alpha0 = internal_scalar(os.path.join(caseDir, "0", "alpha.water"))
        if label == RUNS[0][0] and time == RUNS[0][3]:
            print(f"{'initial':13s} {'0':>5s} {band(alpha0):5d}")
        alpha = internal_scalar(os.path.join(caseDir, time, "alpha.water"))
        air = alpha[(C[:, 2] - Z_INTERFACE) > HEIGHT_ABOVE]
        rhoAir = RHO_AIR + (RHO_WATER - RHO_AIR)*air.mean()
        ratio = RHO_WATER/rhoAir
        print(f"{label:13s} {time:>5s} {band(alpha):5d} {air.size:5d} {int((air > 1e-9).sum()):5d} "
              f"{air.mean():11.3e} {air.max():10.3e} {rhoAir:8.4f} {ratio:7.1f} "
              f"{100*(1 - ratio/(RHO_WATER/RHO_AIR)):11.2f}")


if __name__ == "__main__":
    main()
