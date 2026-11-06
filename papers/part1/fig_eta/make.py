# Copyright (C) 2026 Johan Roenby
# SPDX-License-Identifier: GPL-3.0-or-later
"""Part I, Figure fig:eta: surface elevation (eta_nz80), from ../waveFigures.py.
Output: eta.pdf. Usage: python3 make.py"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
import waveFigures  # noqa: E402

waveFigures.OUT = HERE
waveFigures.fig_eta()
