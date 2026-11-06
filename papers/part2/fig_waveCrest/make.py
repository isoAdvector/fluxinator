# Copyright (C) 2026 Johan Roenby
# SPDX-License-Identifier: GPL-3.0-or-later
"""Part II, Figure fig:waveCrest: crest velocity of the convection variants (waveCrest.pdf).
Made by the shared script ../wave.py, which also writes tab_wave/wave.txt.
Inputs: run/progressiveWave/convection/<case>.
Usage: python3 make.py"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from common import run_script  # noqa: E402

run_script("wave.py")
