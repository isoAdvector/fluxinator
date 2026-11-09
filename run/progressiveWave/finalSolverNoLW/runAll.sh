#!/bin/bash
# Copyright (C) 2026 Johan Roenby
# SPDX-License-Identifier: GPL-3.0-or-later
# LW check: explicit + interfaceUpwind WITHOUT the Lax-Wendroff correction, otherwise identical
# to the fluxinator cases finalSolver_maxCo/{0014,0017,0020,0026,0029,0032} (190x40, 380x80 at
# maxCo 0.5/0.2/0.1, 30 s) and convection/explicitLWIU_190x40_dt0.02 (fixed dt 0.02, diverged with LW).
# Usage: bash runAll.sh   (7 serial runs in parallel)
[ -n "${WM_PROJECT_DIR:-}" ] || { echo "Source the OpenFOAM v2606 environment first" >&2; exit 1; }
cd "$(dirname "$0")"
# The cases 00xx are copies of the fluxinator cases with the correction taken out of the scheme
for c in 0032 0029 0026 0020 0017 0014; do
    [ -d $c ] && continue
    mkdir -p $c/constant
    cp -r ../finalSolver_maxCo/$c/{0.orig,system,Allrun,parmScannerParms} $c/
    for f in $(ls ../finalSolver_maxCo/$c/constant | grep -v polyMesh); do cp -r ../finalSolver_maxCo/$c/constant/$f $c/constant/; done
    sed -i 's/Gauss laxWendroff interfaceUpwind/Gauss interfaceUpwind/' $c/system/fvSchemes
done
for c in 0032 0029 0026 0020 0017 0014 f1_190x40_dt0.02; do
    (cd $c && ./Allrun > log.Allrun 2>&1; echo "$(date +%T) $c $(tail -1 log.fluxinator 2>/dev/null)" >> ../log.runAll) &
done
wait
echo "$(date +%T) all done" >> log.runAll
