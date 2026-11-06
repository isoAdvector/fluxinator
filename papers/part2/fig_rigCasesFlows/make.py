# Copyright (C) 2026 Johan Roenby
# SPDX-License-Identifier: GPL-3.0-or-later
"""Part II, Figure fig:rigCasesFlows. Prescribed-flow cases: ParaView snapshots from snapshots.py
(run that first: pvbatch snapshots.py)
arranged as the figure of the non-uniform flows of Section 5.
Output: rigCasesFlows.pdf. Usage: python3 make.py"""
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

import os
import matplotlib.image as mpimg
from common import plt, save, FULL

SNAP = os.path.join(os.path.dirname(os.path.abspath(__file__)), "snap")


def panels(cases, width, name):
    fig, axes = plt.subplots(1, len(cases), figsize=(width, width/len(cases) + 0.25), squeeze=False)
    axes = axes[0]
    for ax, (case, title) in zip(axes, cases):
        # interpolation "none": the PDF embeds the snapshot at full resolution, not resampled
        ax.imshow(mpimg.imread(os.path.join(SNAP, case + ".png")), extent=(0, 1, 0, 1), interpolation="none")
        ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values():
            s.set_linewidth(0.6); s.set_color("0.4")
        ax.set_title(title, fontsize=9, loc="left")
    fig.tight_layout(w_pad=0.8)
    save(fig, name)


panels([("shear", "(a) shear with cross-flow"), ("rotation", "(b) solid-body rotation")], 0.67*FULL,
       "rigCasesFlows")
