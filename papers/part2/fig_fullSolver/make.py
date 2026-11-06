# Copyright (C) 2026 Johan Roenby
# SPDX-License-Identifier: GPL-3.0-or-later
"""Part II, Figure fig:fullSolver: velocity error of the accelerated disc (fullSolver.pdf).
Made by the shared script ../fullSolver.py, which also writes tab_fullSolver/fullSolver.txt.
Inputs: run/periodicDisc/<case>.
Usage: python3 make.py"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from common import run_script  # noqa: E402

run_script("fullSolver.py")
