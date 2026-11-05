#!/usr/bin/env python3
#
# Writes a summary of each study of the tilted box to data/<study>.txt
# (../common/summariseScan.py: scanned parameters, first p_rgh residual,
# max|U| after the first and last step and its peak), and for meshTimeStep_1step
# also max|U|/dt of the first step.
#
# Usage (from this directory, after the studies have been run):
#   python3 writeData.py

import os
import re
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "data")
STUDIES = ["gravityModels_1step", "gravityModels_100steps", "gravityModels_1000steps",
           "meshTimeStep_1step", "airDensity_1000steps", "nLayers_1000steps"]


def parms(study):
    out = {}
    with open(os.path.join(HERE, study, "parmScannerParms")) as f:
        for line in f:
            m = re.match(r"\S+/(\d{4}):\s+(.*)", line)
            if m:
                out[m.group(1)] = dict(re.findall(r"(\S+): (\S+)", m.group(2)))
    return out


def max_u(study, case):
    f = os.path.join(HERE, study, case, "postProcessing", "fieldMinMax1", "0", "fieldMinMax.dat")
    return np.loadtxt(f, comments="#", usecols=(0, 2), unpack=True, ndmin=2)


def main():
    os.makedirs(OUT, exist_ok=True)
    for study in STUDIES:
        with open(os.path.join(OUT, study + ".txt"), "w") as f:
            subprocess.run([sys.executable, os.path.join(HERE, "..", "common", "summariseScan.py"), study],
                           stdout=f, check=True, cwd=HERE)

    study = "meshTimeStep_1step"
    with open(os.path.join(OUT, study + ".txt"), "a") as f:
        f.write("\nmax|U| after the first time step divided by the time step\n")
        for case, p in parms(study).items():
            t, u = max_u(study, case)
            f.write(f"{case} {p['application']:13s} {p['nx']:>4s}^2 cells  "
                    f"dt={float(p['deltaT']):8.1e}  max|U|={u[0]:.8e} m/s  "
                    f"max|U|/dt={u[0]/float(p['deltaT']):.8e} 1/s\n")
    print("Wrote", ", ".join(sorted(os.listdir(OUT))))


if __name__ == "__main__":
    main()
