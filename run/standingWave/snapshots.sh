#!/bin/bash
#
# Velocity snapshots for the standing-wave figure: interFoam, interIsoFoam and
# interIsoGoam on the refined mesh, written at exactly half and one wave
# period, when the exact velocity field is zero everywhere.
#
# The time step is T/118 (0.010015 s, instead of 0.01 s) so that the writes fall
# exactly on T/2 and T, with T = 2*pi/omega = 1.181816 s from the linear
# dispersion relation for k = pi, d = 0.5 m.
#
# Needs the prepared refined base case of experiments/README.md ($SWBASE,
# default /tmp/swbase).
#
# Usage:  ./snapshots.sh [outdir]      (default: snapshots)
set +u
source ${FOAM_ETC:-/home/roenby/OpenFOAM/OpenFOAM-v2512/etc}/bashrc > /dev/null 2>&1
HERE=$(cd "${0%/*}" && pwd)
OUT=${1:-$HERE/snapshots}
export SWBASE=${SWBASE:-/tmp/swbase}

T=1.181816
DT=$(python3 -c "print($T/118)")

mkdir -p "$OUT"
for app in interFoam interIsoFoam interIsoGoam; do
  $HERE/experiments/runCase.sh "$OUT/$app" app=$app nLayers=2 \
      dt=$DT endTime=$T wint=59 > /dev/null 2>&1
  echo "$app: $(ls -d $OUT/$app/[0-9]* | xargs -n1 basename | tr '\n' ' ')"
done

# Cell centres and volumes for plotting (the mesh is the same for all three)
(cd "$OUT/interIsoGoam" && postProcess -func writeCellCentres -time 0 > /dev/null 2>&1 \
                        && postProcess -func writeCellVolumes -time 0 > /dev/null 2>&1)
