# Copyright (C) 2026 Johan Roenby
# SPDX-License-Identifier: GPL-3.0-or-later
"""Part II, Figure fig:lagWave: time lag in the progressive wave (lagWave.pdf).
Made by the shared script ../lag.py, which also writes fig_lagTGV/lagTGV.pdf and ../data/lag.txt.
Inputs: run/taylorGreenVortex/<case>, ../data/browser_runs.json (browserFiles.py).
Usage: python3 make.py"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from common import run_script  # noqa: E402

run_script("lag.py")
