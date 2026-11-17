#!/bin/bash
# Copyright (C) 2026 Johan Roenby
# SPDX-License-Identifier: GPL-3.0-or-later
# runCT.sh <case> <mode oneStep|carried> <disc origin "x 0 z"> <endTime> "<flow dictionary entries>"
# REF=massCentroid sets the reference point (default cellCentre, see convectionTest.C)
# e.g. runCT.sh rot oneStep "0.5 0 0.75" 1 "Omega 2; rotationCentre (0.5 0 0.5);"
[ -n "${WM_PROJECT_DIR:-}" ] || { echo "Source the OpenFOAM v2606 environment first" >&2; exit 1; }
type dictator > /dev/null 2>&1 || { echo "Source dictator first (github.com/isoAdvector/dictator)" >&2; exit 1; }
out=$1; mode=$2; origin=$3; tEnd=$4; flow=$5
BASE=$(readlink -f "$(dirname "$0")/baseCase")
rm -rf "$out"; mkdir -p "$out"; cp -r $BASE/* "$out"/; cd "$out" || exit 1
set -e   # stop if dictator fails (a missing file or an ambiguous entry)
{
[ -n "$N" ] && dictator system/blockMeshDict n $N
[ -n "$DT" ] && dictator system/controlDict deltaT $DT
dictator system/setAlphaFieldDict origin "($origin)"
dictator system/controlDict endTime $tEnd
} > log.dictator
set +e
cat > system/convectionTestDict <<EOD
FoamFile { version 2.0; format ascii; class dictionary; object convectionTestDict; }
mode $mode;
${REF:+reference $REF;}
flow { $flow }
#include "variants"
EOD
./Allrun > log.Allrun 2>&1
grep -q "^End" log.convectionTest || echo "FAILED $out"
