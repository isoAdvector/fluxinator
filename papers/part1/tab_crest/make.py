"""Part I, Table tab:crest: crest kinematics (../crestKinematics.py). Output: tab_crest.txt. Usage: python3 make.py"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
import crestKinematics  # noqa: E402

crestKinematics.table(os.path.join(HERE, "tab_crest.txt"))
