# Copyright (C) 2026 Johan Roenby
# SPDX-License-Identifier: GPL-3.0-or-later
"""Time-lag toy: translating Taylor-Green vortex, single phase, doubly periodic. The exact
x-mean <u>(z) stays U0. With the mass flux at t^n and the convected velocity at t^n + theta dt
the predicted spurious mean flow is <u>(z, t) - U0 = a(t) cos 2z, a(t) = -theta dt U0 t / 2.
Per case: a(t) = (2/N) sum_z (<u> - U0) cos 2z, the total mean <u> - U0 (momentum check) and
the prediction for the case's theta (file theta in the case). Writes tgv.txt and tgv.png (not used by the paper figures).
Usage: python3 writeData.py <case> [<case> ...]"""
import os
import re
import sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "common"))
import foamio  # noqa: E402
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
n = 64
zc = (np.arange(n) + 0.5)*2*np.pi/n
fig, ax = plt.subplots(figsize=(5.0, 3.4))
L = [f"{'case':26s} {'theta':>5s} {'dt':>7s} {'a(10) meas':>11s} {'a(10) pred':>11s} {'ratio':>6s} {'<u>-U0':>9s}"]
for case in sys.argv[1:]:
    U0 = float(open(os.path.join(case, "U0")).read())
    th = float(open(os.path.join(case, "theta")).read())
    dt = float(foamio.entry(os.path.join(case, "system/controlDict"), "deltaT"))
    ts = sorted((d for d in os.listdir(case) if re.fullmatch(r"\d+(\.\d+)?", d)), key=float)
    t, a, m = [], [], []
    for d in ts:
        U = foamio.read_vol(os.path.join(case, d, "U")).reshape(n, n, 3)
        um = U[:, :, 0].mean(1) - U0
        t.append(float(d)); a.append(2*np.mean(um*np.cos(2*zc))); m.append(U[:, :, 0].mean() - U0)
    pred = -0.5*th*dt*U0*t[-1]
    line, = ax.plot(t, a, marker="o", ms=3, label=os.path.basename(case.rstrip("/")))
    ax.plot(t, [-0.5*th*dt*U0*x for x in t], "--", color=line.get_color(), lw=0.8)
    L.append(f"{os.path.basename(case.rstrip('/')):26s} {th:5.2f} {dt:7.4f} {a[-1]:+11.5f} {pred:+11.5f} "
             f"{(a[-1]/pred if pred else float('nan')):6.2f} {m[-1]:+9.2e}")
ax.set(xlabel="t", ylabel=r"amplitude of $\cos 2z$ in $\langle u\rangle - U_0$")
ax.axhline(0, color="0.4", lw=0.7)
ax.legend(fontsize=7)
ax.set_title("solid: measured; dashed: prediction", fontsize=9)
fig.tight_layout()
fig.savefig(os.path.join(HERE, "tgv.png"), dpi=150)
open(os.path.join(HERE, "tgv.txt"), "w").write("\n".join(L) + "\n")
print("\n".join(L))
