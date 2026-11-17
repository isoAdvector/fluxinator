#!/bin/bash
# Copyright (C) 2026 Johan Roenby
# SPDX-License-Identifier: GPL-3.0-or-later
# Where the largest velocity of the standing wave is: the two-period runs of paperScan.sh
# (plicRDF) with the fields written every time step, split by airWater.py into air
# (alpha < 0.01), interface cells and water (alpha > 0.99). Writes <out>/airWater.csv.
# Usage: ./airWater.sh [<out>]   (default solvers_airWater)
set +u
[ -n "${WM_PROJECT_DIR:-}" ] || { echo "Source the OpenFOAM v2606 environment first" >&2; exit 1; }
HERE=$(cd "${0%/*}" && pwd)
OUT=${1:-$HERE/solvers_airWater}

mkdir -p "$OUT"
for mesh in uniform refined; do
  case $mesh in
    uniform) base=$HERE/baseUniformMesh ;;
    refined) base=$HERE/baseRefinedMesh ;;
  esac
  for app in interFoam interIsoFoam interIsoGoam; do
    SWBASE=$base $HERE/runCase.sh "$OUT/${mesh}_$app" app=$app endTime=2.4 wint=1 nLayers=2 \
        recon=plicRDF > /dev/null 2>&1 &
  done
done
wait
python3 "$HERE/airWater.py" "$OUT"
