"""Part II, Figure fig:selCost: cost of the selected runs (selCost.pdf).
Made by the shared script ../selection.py, which also writes fig_selCase, fig_selHeight and ../data/paperSelection.txt.
Inputs: ../data/browser.json (python3 ../browserData.py) and run/progressiveWave/data/timing_maxCo0.1.csv.
Usage: python3 make.py"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from common import run_script  # noqa: E402

run_script("selection.py")
