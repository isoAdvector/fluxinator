#!/bin/bash
# Copyright (C) 2026 Johan Roenby
# SPDX-License-Identifier: GPL-3.0-or-later
# runStored.sh <storedCase> <outdir> [nProcs]
# Rerun a stored Part I case (run/<case>/<study>/<run>, any Allrun-driven case)
# with the single change PIMPLE ddtCorr no (DDTCORR=yes in the environment for
# regression checks). Stock solvers are
# replaced by copies with the ddtCorr switch (applications/test/ddtCorrSwitch):
# interFoam -> interFoamDdtCorr, interIsoFoam -> interIsoFoamDdtCorr; with
# ddtCorr yes these reproduce the stored runs bitwise (reg/). Figure rendering
# (pvpython) is skipped. With nProcs > 1 the wave case runs with AllrunParallel
# on nProcs subdomains and is reconstructed afterwards.
[ -n "${WM_PROJECT_DIR:-}" ] || { echo "Source the OpenFOAM v2606 environment first" >&2; exit 1; }
type dictator > /dev/null 2>&1 || { echo "Source dictator first (github.com/isoAdvector/dictator)" >&2; exit 1; }
here=$(dirname "$(readlink -f "$0")")   # before the cd below, which breaks a relative $0
S=$(readlink -f "$1"); out=$2; np=${3:-1}
rm -rf "$out"; mkdir -p "$out"
for f in $S/*; do
    b=$(basename $f)
    case $b in
        processor*|postProcessing|log.*|*.png|*.foam|dynamicCode) continue;;
    esac
    [[ $b =~ ^[0-9.e+-]+$ ]] && continue     # time directories
    cp -r $f "$out"/
done
rm -rf "$out"/constant/polyMesh
cd "$out" || exit 1
sed -i 's/^\( *\)pvpython/\1: pvpython/' Allrun
set -e   # stop if dictator fails (a missing file or an ambiguous entry)
{
for d in system/controlDict system/controlDict.run; do
    [ -f $d ] || continue
    app=$(dictator $d application)
    case $app in
        interFoam)    dictator $d application interFoamDdtCorr;;
        interIsoFoam) dictator $d application interIsoFoamDdtCorr;;
    esac
done
dictator system/fvSolution PIMPLE.ddtCorr ${DDTCORR:-no}
} > log.dictator
set +e
echo "stored: $1  change: PIMPLE.ddtCorr ${DDTCORR:-no}  nProcs: $np" > ddtCorrOff.info
if [ "$np" -gt 1 ]; then
    cp "$here"/../progressiveWave/baseCase_fixedDt/AllrunParallel .
    set -e
    {
    dictator system/decomposeParDict numberOfSubdomains $np
    dictator system/decomposeParDict method simple
    dictator system/decomposeParDict simpleCoeffs.n "($np 1 1)"
    } >> log.dictator
    set +e
    ./AllrunParallel > log.Allrun 2>&1
    # The fixed-dt cases write no fields, so there is nothing to reconstruct
    { reconstructPar > log.reconstructPar 2>&1 || grep -q 'No times selected' log.reconstructPar; } \
        && rm -rf processor*
else
    ./Allrun > log.Allrun 2>&1
fi
app=$(dictator system/controlDict application)
grep -q "^End" log.$app || echo "FAILED $out"
