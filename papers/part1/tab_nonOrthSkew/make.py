#!/usr/bin/env python3
"""Part I, Table tab:nonOrthSkew: maximum velocity magnitude [m/s] in the skewed box
after 10 time steps (t = 0.1 s). Each entry names the run it comes from
(run/skewedBox/<study>/<case>). Output: tab_nonOrthSkew.txt. Usage: python3 make.py"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
from common import case, max_u_at, write_table  # noqa: E402

COLS = ["corrected, n_c = 0", "corrected, n_c = 2", "corrected, hydrostatic",
        "uncorrected", "uncorrected, rows level"]
# row -> (study, case) per column; None: not run
ROWS = {
    "interFoam, basic": [("schemes_10steps", "0000"), ("nonOrthCorrectors_10steps", "0009"),
                         ("hydrostaticStart_10steps", "0005"), ("schemes_10steps", "0003"),
                         ("schemesRotatedMesh_10steps", "0003")],
    "interIsoFoam, basic": [("schemes_10steps", "0001"), ("nonOrthCorrectors_10steps", "0010"),
                            ("hydrostaticStart_10steps", "0004"), ("schemes_10steps", "0004"),
                            ("schemesRotatedMesh_10steps", "0004")],
    "interIsoGoam, constantHeight": [("nonOrthCorrectors_10steps", "0007"), ("nonOrthCorrectors_10steps", "0008"),
                                     ("hydrostaticStart_10steps", "0002"), ("schemes_10steps", "0006"), None],
    "interIsoGoam, geometric, nL = 1": [("schemes_10steps", "0002"), ("nonOrthCorrectors_10steps", "0002"),
                                        ("hydrostaticStart_10steps", "0000"), ("schemes_10steps", "0005"),
                                        ("schemesRotatedMesh_10steps", "0005")],
}

out = ["Maximum velocity magnitude [m/s] in the skewed box after 10 time steps (t = 0.1 s)",
       f"{'':32s}" + "".join(f" {c:>24s}" for c in COLS)]
src = ["", "Sources (run/skewedBox/<study>/<case>):"]
for row, cells in ROWS.items():
    vals = []
    for col, cell in zip(COLS, cells):
        if cell is None:
            vals.append(f" {'--':>24s}")
            continue
        study, c = cell
        vals.append(f" {max_u_at(case('skewedBox', study), c, 0.1):24.2e}")
        src.append(f"  {row}, {col}: {study}/{c}")
    out.append(f"{row:32s}" + "".join(vals))
write_table(HERE, "tab_nonOrthSkew", out + src)
