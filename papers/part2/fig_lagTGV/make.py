"""Part II, Figure fig:lagTGV: time lag in the translating Taylor-Green vortex (lagTGV.pdf).
Made by the shared script ../lag.py, which also writes fig_lagWave/lagWave.pdf and ../data/lag.txt.
Inputs: run/taylorGreenVortex/<case>, ../data/browser_runs.json (browserFiles.py).
Usage: python3 make.py"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from common import run_script  # noqa: E402

run_script("lag.py")
