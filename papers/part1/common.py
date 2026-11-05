"""Shared paths for the Part I figure and table scripts.

Every figure or table of Part I with LaTeX label fig:X or tab:X has its own
directory fig_X/ or tab_X/ here, with a make.py (or make.sh) that writes the
file the paper shows into that directory.
"""
import os
import sys

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
RUN = os.path.join(REPO, "run")
sys.path.insert(0, os.path.join(RUN, "common"))


def case(name, *rest):
    """Path in run/<name>/..."""
    return os.path.join(RUN, name, *rest)


def max_u(study_dir, case_name):
    """Time and max|U| of a case from its fieldMinMax1 function object."""
    import numpy as np
    f = os.path.join(study_dir, case_name, "postProcessing", "fieldMinMax1", "0", "fieldMinMax.dat")
    return np.loadtxt(f, comments="#", usecols=(0, 2), unpack=True, ndmin=2)


def max_u_at(study_dir, case_name, t):
    """max|U| of a case at the output time nearest to t."""
    import numpy as np
    tt, u = max_u(study_dir, case_name)
    return u[np.argmin(abs(tt - t))]


def write_table(here, name, lines):
    """Write the table file <name>.txt in the table directory and print it."""
    open(os.path.join(here, name + ".txt"), "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))
