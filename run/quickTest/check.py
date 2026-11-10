# Copyright (C) 2026 Johan Roenby
# SPDX-License-Identifier: GPL-3.0-or-later
"""Prints the results of the quick test next to the expected values."""
import glob
import os
import re

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def latest(case):
    times = [d for d in glob.glob(os.path.join(HERE, case, "[0-9]*")) if re.fullmatch(r"[0-9.e+-]+", os.path.basename(d))]
    return max(times, key=lambda d: float(os.path.basename(d)))


def internal(path, kind):
    """The internal field of a volScalarField or volVectorField file in ascii format."""
    m = re.search(r"internalField\s+nonuniform List<%s>\s*\d+\s*\((.*?)\n\)\s*;" % kind, open(path).read(), re.S)
    if kind == "scalar":
        return np.array(m.group(1).split(), dtype=float)
    return np.array([v.split() for v in re.findall(r"\(([^()]*)\)", m.group(1))], dtype=float)


def finished(case):
    log = os.path.join(HERE, case, "log.fluxinator")
    return os.path.exists(log) and re.search(r"^End", open(log).read(), re.M) is not None


if finished("disc"):
    t = latest("disc")
    U = internal(os.path.join(t, "U"), "vector")
    exact = np.array([0.1 + 0.5*float(os.path.basename(t)), 0, 0.05])
    print(f"disc: t = {os.path.basename(t)} s, max|U - U0 - a t| = {np.abs(U - exact).max():.1e} m/s (expected below 1e-8)")
else:
    # Before OpenFOAM v2506 plicRDF stops on the cyclic patches of the disc
    # (bin/adjustCodeToLoadedOFVersion)
    need = [l.strip() for f in glob.glob(os.path.join(HERE, "disc", "log.*")) for l in open(f, errors="replace")
            if "needs OpenFOAM v2506" in l]
    print(f"disc: stopped: {need[0]}" if need else "disc: did not finish, see disc/log.*")
if os.path.isdir(os.path.join(HERE, "wave")):
    if finished("wave"):
        t = latest("wave")
        u = internal(os.path.join(t, "U"), "vector")[:, 0]
        water = internal(os.path.join(t, "alpha.water"), "scalar") > 0.5
        print(f"wave: t = {os.path.basename(t)} s, largest horizontal water velocity in a cell = {u[water].max():.3f} m/s "
              "(expected 0.426 m/s; exact crest velocity 0.445 m/s)")
    else:
        print("wave: did not finish, see wave/log.*")
