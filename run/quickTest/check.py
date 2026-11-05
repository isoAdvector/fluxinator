"""Prints the results of the quick test next to the expected values."""
import glob
import os
import re

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def latest(case):
    return max((d for d in glob.glob(os.path.join(HERE, case, "[0-9]*")) if os.path.isdir(d)),
               key=lambda d: float(os.path.basename(d)))


def vectors(path):
    """The internal field of a volVectorField file in ascii format."""
    m = re.search(r"internalField\s+nonuniform List<vector>\s*\d+\s*\((.*?)\n\)\s*;", open(path).read(), re.S)
    return np.array([v.split() for v in re.findall(r"\(([^()]*)\)", m.group(1))], dtype=float)


def finished(case):
    log = os.path.join(HERE, case, "log.interIsoGoam")
    return os.path.exists(log) and re.search(r"^End", open(log).read(), re.M) is not None


if finished("disc"):
    t = latest("disc")
    U = vectors(os.path.join(t, "U"))
    exact = np.array([0.1 + 0.5*float(os.path.basename(t)), 0, 0.05])
    print(f"disc: t = {os.path.basename(t)} s, max|U - U0 - a t| = {np.abs(U - exact).max():.1e} m/s (expected below 1e-8)")
else:
    print("disc: did not finish, see disc/log.*")
if os.path.isdir(os.path.join(HERE, "wave")):
    if finished("wave"):
        t = latest("wave")
        print(f"wave: t = {os.path.basename(t)} s, max|U| = {np.linalg.norm(vectors(os.path.join(t, 'U')), axis=1).max():.3f} m/s "
              "(expected EXPECTEDWAVE; exact crest speed 0.445 m/s)")
    else:
        print("wave: did not finish, see wave/log.*")
