"""Part I, Figure fig:spurious: maximum velocity (spurious_velocity), from ../waveFigures.py.
Output: spurious.pdf. Usage: python3 make.py"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
import waveFigures  # noqa: E402

waveFigures.OUT = HERE
waveFigures.fig_umax()
