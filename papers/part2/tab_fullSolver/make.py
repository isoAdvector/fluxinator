"""Part II, Table tab:fullSolver: accelerated disc and uniform flow, full solver (fullSolver.txt).
Made by the shared script ../fullSolver.py, which also writes fig_fullSolver/fullSolver.pdf.
Inputs: run/periodicDisc/<case>.
Usage: python3 make.py"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from common import run_script  # noqa: E402

run_script("fullSolver.py")
