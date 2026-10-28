#!/usr/bin/env python3
#
# Writes the numbers behind the tilted box results in the paper to paperData/
# (tracked in git, unlike the scans themselves):
#   maxU_gravityModelScan.txt  max|U| at t = 0.01, 0.1, 1, 10 s and its peak
#                              (gravityModelScan_1000deltaT)
#   meshDeltaTScan.txt         max|U| after the first step and max|U|/dt, for
#                              four meshes, five time steps and both solvers
#   <scan>.txt                 ../summariseScan.py output for each scan
#   interfaceBand.txt          analyseInterfaceBand.py output
#   vertexOffsetTheory.txt      closed-form prediction of which cells are hit,
#                              checked against the interFoam runs of
#                              tiltScan_1deltaT and waterLevelScan_1deltaT
#
# Usage (from this directory, after the scans have been run):
#   python3 paperData.py

import os
import re
import subprocess
import sys

import numpy as np

OUT = "paperData"
SCANS = ["gravityModelScan_1deltaT", "gravityModelScan_100deltaT",
         "gravityModelScan_1000deltaT", "meshDeltaTScan",
         "airDensityScan_1000deltaT",
         "nLayersScan", "tiltScan_1deltaT", "waterLevelScan_1deltaT"]


def parms(scan):
    out = {}
    with open(os.path.join(scan, "parmScannerParms")) as f:
        for line in f:
            m = re.match(r"\S+/(\d{4}):\s+(.*)", line)
            if m:
                out[m.group(1)] = dict(re.findall(r"(\S+): (\S+)", m.group(2)))
    return out


def max_u(scan, case):
    f = os.path.join(scan, case, "postProcessing", "fieldMinMax1", "0", "fieldMinMax.dat")
    return np.loadtxt(f, comments="#", usecols=(0, 2), unpack=True, ndmin=2)


def run(cmd, outFile):
    with open(os.path.join(OUT, outFile), "w") as f:
        subprocess.run(cmd, stdout=f, check=True)


def main():
    os.makedirs(OUT, exist_ok=True)

    for scan in SCANS:
        run([sys.executable, "../summariseScan.py", scan], scan + ".txt")
    run([sys.executable, "analyseInterfaceBand.py"], "interfaceBand.txt")
    run([sys.executable, "vertexOffsetTheory.py", "--validate",
         "tiltScan_1deltaT", "waterLevelScan_1deltaT"], "vertexOffsetTheory.txt")

    scan = "gravityModelScan_1000deltaT"
    with open(os.path.join(OUT, "maxU_gravityModelScan.txt"), "w") as f:
        f.write(f"max|U| [m/s] in {scan}\n")
        f.write(f"{'case':5s} {'gravityModel':15s} {'application':13s}"
                + "".join(f"{'t=' + str(T):>11s}" for T in (0.01, 0.1, 1, 10)) + f"{'peak':>11s}\n")
        for case, p in parms(scan).items():
            t, u = max_u(scan, case)
            vals = [u[np.argmin(abs(t - T))] for T in (0.01, 0.1, 1, 10)]
            f.write(f"{case:5s} {p['gravityModel']:15s} {p['application']:13s}"
                    + "".join(f"{v:11.3e}" for v in vals) + f"{u.max():11.3e}\n")

    scan = "meshDeltaTScan"
    with open(os.path.join(OUT, "meshDeltaTScan.txt"), "a") as f:
        f.write("\nmax|U| after the first time step divided by the time step\n")
        for case, p in parms(scan).items():
            t, u = max_u(scan, case)
            f.write(f"{case} {p['application']:13s} {p['nx']:>4s}^2 cells  "
                    f"dt={float(p['deltaT']):8.1e}  max|U|={u[0]:.8e} m/s  "
                    f"max|U|/dt={u[0]/float(p['deltaT']):.8e} 1/s\n")
    print("Wrote", ", ".join(sorted(os.listdir(OUT))))


if __name__ == "__main__":
    main()
