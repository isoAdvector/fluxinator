#!/bin/bash
#
# Figure for the standing-wave section: ParaView snapshots of the three solvers
# at t = T (run ./snapshots.sh first), assembled by figure_standingWave.tex.
#
# Usage:  ./figure.sh [time] [scale]     (defaults: 1.181817, 0.1 m per m/s)
set -e
cd "${0%/*}"
TIME=${1:-1.181817}
SCALE=${2:-0.1}
for app in interFoam interIsoFoam interIsoGoam; do
  ref=""; [ $app = interFoam ] && ref="--reference"
  pvbatch generateSnapshotFigure.py snapshots/$app $TIME snapshots/$app.png \
          --scale $SCALE $ref > snapshots/log.$app.png 2>&1
done
pdflatex -interaction=nonstopmode figure_standingWave.tex > log.figure_standingWave 2>&1
echo "wrote figure_standingWave.pdf"
