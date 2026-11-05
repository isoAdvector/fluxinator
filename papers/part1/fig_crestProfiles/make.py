"""Part I, Figure fig:crestProfiles: velocity under the crest (../crestKinematics.py). Output: crestProfiles.pdf. Usage: python3 make.py"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
import crestKinematics  # noqa: E402

crestKinematics.profiles(os.path.join(HERE, "crestProfiles"))
