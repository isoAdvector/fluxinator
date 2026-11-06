# Copyright (C) 2026 Johan Roenby
# SPDX-License-Identifier: GPL-3.0-or-later
"""Part II, Table tab:rigConvergence: mesh convergence of the prescribed flows (rigConvergence.txt).
Made by the shared script ../rig.py, which also writes ../data/rig.txt.
Inputs: run/prescribedFlows/<case>_n{25,100}.
Usage: python3 make.py"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from common import run_script  # noqa: E402

run_script("rig.py")
