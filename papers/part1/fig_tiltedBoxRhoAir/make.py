#!/usr/bin/env python
#
# Plots the maximum velocity magnitude at the end of the runs against the air
# density for interFoam, interIsoFoam (interIsoGoam with the basic gravity
# model) and interIsoGoam with the geometric gravity model in the tilted box.
# Part I, Figure fig:tiltedBoxRhoAir. Input: run/tiltedBox/airDensity_1000steps
# (t = 10 s). Runs that stopped before the end time are drawn with a cross at
# their last value. Output: tiltedBoxRhoAir.pdf (and .png).
#
# Usage:
#   python3 make.py

import argparse
import os
import re
import sys

import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
from common import case  # noqa: E402

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 10,
    "axes.linewidth": 0.8,
    "lines.linewidth": 1.4,
})

# (application, gravityModel) -> label, colour, marker
STYLE = {
    ("interFoam", "basic"): ("interFoam, basic", "#2a78d6", "o"),
    ("interIsoGoam", "basic"): ("interIsoFoam, basic", "#eb6834", "s"),
    ("interIsoGoam", "geometric"): ("interIsoGoam, geometric", "#1baf7a", "D"),
}
RHO_WATER = 1000.0
Y_RANGE = (1e-12, 1e6)

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--scanDir", default=case("tiltedBox", "airDensity_1000steps"))
parser.add_argument("--endTime", type=float, default=10.0, help="End time of the scan runs [s]")
parser.add_argument("--out", default=os.path.join(HERE, "tiltedBoxRhoAir"), help="Output file name without extension")
args = parser.parse_args()
END_TIME = args.endTime

data = {}
with open(os.path.join(args.scanDir, "parmScannerParms")) as f:
    for line in f:
        case = line.split(":")[0].split("/")[-1]
        parms = dict(re.findall(r"([\w.]+):\s*(\S+)", line.split(":", 1)[1]))
        key = (parms["application"], parms["gravityModel"])
        fmm = os.path.join(args.scanDir, case, "postProcessing", "fieldMinMax1", "0", "fieldMinMax.dat")
        if not os.path.isfile(fmm):
            print(f"{case}: no fieldMinMax output")
            continue
        t, umax = np.loadtxt(fmm, comments="#", unpack=True, usecols=(0, 2), ndmin=2)
        complete = t[-1] >= END_TIME - 1e-9 and np.isfinite(umax[-1])
        rho = float(parms["air.rho"])
        if rho >= RHO_WATER:
            # No density jump: velocities are at round-off for every model
            print(f"{case}: rho_air = {rho:g} (no density jump) max|U| = {umax[-1]:.3e}, not plotted")
            continue
        data.setdefault(key, []).append((rho, float(umax[-1]), bool(complete), case, t[-1]))

fig, ax = plt.subplots(figsize=(4.6, 3.4))
for key, (label, colour, marker) in STYLE.items():
    pts = sorted(data.get(key, []))
    print(label)
    for rho, u, ok, case, tl in pts:
        print(f"  {case} rho_air={rho:8.1e}  max|U|={u:.3e}  {'t=1' if ok else f'STOPPED at t={tl}'}")
    ok_pts = [(r, u) for r, u, ok, *_ in pts if ok]
    if ok_pts:
        ax.loglog([p[0] for p in ok_pts], [p[1] for p in ok_pts], "-", color=colour,
                  marker=marker, markersize=6, label=label)
    bad = [(r, u) for r, u, ok, *_ in pts if not ok]
    if bad:
        # Diverged runs: cross kept just inside the top of the axis
        ax.loglog([p[0] for p in bad], [min(max(p[1], Y_RANGE[0]), 0.3*Y_RANGE[1]) for p in bad], "x",
                  color=colour, markersize=8, markeredgewidth=1.8)
ax.set_xlabel(r"$\rho_{air}$ [kg/m$^3$] ($\rho_{water} = 1000$ kg/m$^3$)")
ax.set_ylabel(r"max$|\mathbf{U}|$ at $t = %g$ s [m/s]" % END_TIME)
ax.set_ylim(*Y_RANGE)
ax.grid(True, which="major", linewidth=0.4, alpha=0.6)
ax.legend(loc="upper right", fontsize=8, frameon=False)
fig.tight_layout()
fig.savefig(args.out + ".pdf")
fig.savefig(args.out + ".png", dpi=200)
print("Wrote", args.out + ".pdf and .png")
