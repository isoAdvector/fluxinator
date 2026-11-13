#!/bin/bash
# Copyright (C) 2026 Johan Roenby
# SPDX-License-Identifier: GPL-3.0-or-later
#
# Velocity snapshots for the standing-wave figure: interFoam, interIsoFoam and
# interIsoGoam on the refined mesh, written at exactly half and one wave
# period, when the exact velocity field is zero everywhere.
#
# The time step is T/118 (0.010015 s, instead of 0.01 s) so that the writes fall
# exactly on T/2 and T, with T = 2*pi/omega = 1.181816 s from the linear
# dispersion relation for k = pi, d = 0.5 m.
#
# Uses the prepared refined base case baseRefinedMesh (mesh built; $SWBASE overrides).
#
# Usage:  ./snapshots.sh [outdir]      (default: snapshots)
#         RECON=plicRDF ./snapshots.sh snapshots_plicRDF   (reconstruction, default isoAlpha)
set +u
[ -n "${WM_PROJECT_DIR:-}" ] || { echo "Source the OpenFOAM v2606 environment first" >&2; exit 1; }
HERE=$(cd "${0%/*}" && pwd)
OUT=${1:-$HERE/snapshots}
export SWBASE=${SWBASE:-$HERE/baseRefinedMesh}
RECON=${RECON:-isoAlpha}

T=1.181816
DT=$(python3 -c "print($T/118)")

mkdir -p "$OUT"
for app in interFoam interIsoFoam interIsoGoam; do
  $HERE/runCase.sh "$OUT/$app" app=$app nLayers=2 recon=$RECON \
      dt=$DT endTime=$T wint=59 > /dev/null 2>&1
  echo "$app: $(ls -d $OUT/$app/[0-9]* | xargs -n1 basename | tr '\n' ' ')"
done

# Cell centres and volumes for plotting (the mesh is the same for all three)
(cd "$OUT/interIsoGoam" && postProcess -func writeCellCentres -time 0 > /dev/null 2>&1 \
                        && postProcess -func writeCellVolumes -time 0 > /dev/null 2>&1)
