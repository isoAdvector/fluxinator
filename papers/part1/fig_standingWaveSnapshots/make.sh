#!/bin/bash
# Copyright (C) 2026 Johan Roenby
# SPDX-License-Identifier: GPL-3.0-or-later
# Part I, Figure fig:standingWaveSnapshots: velocity snapshots of the standing wave at
# t = T for interFoam, interIsoFoam and interIsoGoam.
# Inputs: run/standingWave/snapshots/<solver>/1.181817 (run/standingWave/snapshots.sh).
# Output: standingWaveSnapshots.pdf (and the panels <solver>.png).
#
# Usage:  ./make.sh [time] [scale]     (defaults: 1.181817, 0.1 m per m/s)
set -e
cd "${0%/*}"
TIME=${1:-1.181817}
SCALE=${2:-0.1}
SNAP=../../../run/standingWave/snapshots
for app in interFoam interIsoFoam interIsoGoam; do
  ref=""; [ $app = interFoam ] && ref="--reference"
  pvbatch snapshotImage.py $SNAP/$app $TIME $app.png \
          --scale $SCALE $ref > log.$app 2>&1
done
pdflatex -interaction=nonstopmode standingWaveSnapshots.tex > log.make 2>&1
rm -f standingWaveSnapshots.aux standingWaveSnapshots.log
