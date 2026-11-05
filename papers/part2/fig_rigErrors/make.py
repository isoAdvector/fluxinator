"""Kinematic rig: largest velocity error per time step against time for the convection
variants. (a) shear with cross-flow, one-step mode (error of one step from the exact field);
(b) solid-body rotation, carried mode (error accumulated over the run). 50x50, dt 0.005 s.
Part II, Figure fig:rigErrors. Input: run/prescribedFlows/{shear,rotationCarried}. Output: rigErrors.pdf.
Usage: python3 make.py"""
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

import os
import numpy as np
from common import RIG, plt, save, FULL, COL

VAR = [("implEuler_new", "implicit", COL["implicit"], "-"),
       ("explEuler_new", "explicit", COL["explicit"], "-"),
       ("explIU_new", "explicit + interfaceUpwind", COL["explicit + interfaceUpwind"], "-")]
fig, axes = plt.subplots(1, 2, figsize=(FULL, 2.7))
for ax, (case, title) in zip(axes, [("shear", "(a) shear with cross-flow, one step"),
                                     ("rotationCarried", "(b) solid-body rotation, carried")]):
    for v, lab, c, ls in VAR:
        d = np.loadtxt(os.path.join(RIG, case, "postProcessing", "convectionTest", v + ".dat"), ndmin=2)
        ax.semilogy(d[:, 0], d[:, 1], ls, color=c, lw=1.0, label=lab)
    ax.set_title(title, fontsize=9, loc="left"); ax.set_xlabel("t [s]")
axes[0].set_ylabel("max |U - u| [m/s]")
axes[1].legend(loc="center right", fontsize=8, framealpha=0.9)
fig.tight_layout()
save(fig, "rigErrors")
