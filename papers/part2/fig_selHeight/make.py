"""Part II, Figure fig:selHeight: wave height of the selected runs (selHeight.pdf).
Made by the shared script ../selection.py, which also writes fig_selCase, fig_selCost and ../data/paperSelection.txt.
Inputs: ../data/browser.json (python3 ../browserData.py).
Usage: python3 make.py"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from common import run_script  # noqa: E402

run_script("selection.py")
