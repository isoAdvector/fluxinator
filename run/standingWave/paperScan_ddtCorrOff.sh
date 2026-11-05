#!/bin/bash
#
# Copy of paperScan.sh with PIMPLE ddtCorr no (runCase_ddtCorrOff.sh).
# Data for the standing-wave section: spurious air velocities of interFoam,
# interIsoFoam and interIsoGoam on a uniform and on a locally refined mesh.
#
# The wave starts at rest at maximum elevation, so the exact velocity after one
# step of 0.01 s is only 0.0154 m/s and anything above that is error. The error
# depends strongly on where the interface falls relative to the cell faces, so
# every entry is repeated over NOFF interface offsets spanning one coarse cell
# and the median is reported.
#
# Usage:  ./paperScan.sh [outdir]      (default: scan)
#
# Writes <outdir>/scan.csv and <outdir>/history_<app>.csv
set +u
[ -n "${WM_PROJECT_DIR:-}" ] || { echo "Source the OpenFOAM v2512 environment first" >&2; exit 1; }
HERE=$(cd "${0%/*}" && pwd)
OUT=${1:-$HERE/solvers_ddtCorrOff}
NOFF=${NOFF:-12}
DZ=0.02                      # coarse cell height
RUN=$HERE/runCase_ddtCorrOff.sh

mkdir -p "$OUT"
csv=$OUT/scan.csv
echo "mesh,app,offset,first,peak20" > $csv

for mesh in uniform refined; do
  case $mesh in
    uniform) export SWBASE=$HERE/baseUniformMesh ;;
    refined) export SWBASE=$HERE/baseRefinedMesh ;;
  esac
  for app in interFoam interIsoFoam interIsoGoam; do
    for i in $(seq 0 $((NOFF-1))); do
      z=$(python3 -c "print(0.5 + $i*$DZ/$NOFF)")
      d=$OUT/${mesh}_${app}_$i
      $RUN "$d" app=$app endTime=0.2 wint=20 nLayers=2 \
           alphaOrigin="(-0.5 0 $z)" > /dev/null 2>&1
      first=$(grep "min/max(mag(U))" $d/log.run | head -1 | awk '{print $NF}')
      peak=$(grep "min/max(mag(U))" $d/log.run | awk '{print $NF}' | sort -g | tail -1)
      echo "$mesh,$app,$z,$first,$peak" >> $csv
    done
  done
done

# Time history over two wave periods on the refined mesh, middle offset
export SWBASE=$HERE/baseRefinedMesh
for app in interFoam interIsoFoam interIsoGoam; do
  d=$OUT/hist_$app
  $RUN "$d" app=$app endTime=2.4 wint=120 nLayers=2 > /dev/null 2>&1
  awk '/^Time = /{t=$3} /min\/max\(mag\(U\)\)/{print t","$NF}' $d/log.run \
      > $OUT/history_$app.csv
done
echo "wrote $csv and $OUT/history_*.csv"
