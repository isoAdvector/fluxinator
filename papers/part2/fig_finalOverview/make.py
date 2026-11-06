# Copyright (C) 2026 Johan Roenby
# SPDX-License-Identifier: GPL-3.0-or-later
"""Part II, Figure fig:finalOverview: fluxinator against the Part I solvers rerun with ddtCorr
off (finalOverview.pdf), with the numbers in finalSolver.txt and finalSolver.csv.
Made by ../finalSolver.py (main). Inputs: run/progressiveWave/{finalSolver,solvers}_{maxCo,fixedDt}*.
Usage: python3 make.py [-j N]"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import finalSolver  # noqa: E402

if __name__ == "__main__":
    finalSolver.main(int(sys.argv[sys.argv.index("-j") + 1]) if "-j" in sys.argv else 8)
