#!/bin/bash
# Copyright (C) 2026 Johan Roenby
# SPDX-License-Identifier: GPL-3.0-or-later
# Builds the timing cases: three repeats per solver and mesh, from the finished set-ups of
# solvers_maxCo_ddtCorrOff and finalSolver_maxCo, with endTime 6 and the wave gauge off.
cd "${0%/*}" || exit 1
type dictator > /dev/null 2>&1 || { echo "Source dictator first (github.com/isoAdvector/dictator)" >&2; exit 1; }
wave=..
make1() {   # make1 <source case> <name>
    for r in 1 2 3; do
        d=$2_r$r
        mkdir -p "$d/constant"
        cp -r "$1"/{0.orig,system,Allrun} "$d"/
        for f in $(ls "$1/constant" | grep -v polyMesh); do cp -r "$1/constant/$f" "$d/constant/"; done
        {
        dictator "$d/system/controlDict.run" endTime 6 &&
        dictator "$d/system/controlDict.run" functions.interfaceHeight1.enabled false
        } > "$d/log.dictator" || exit 1
    done
}
off=$wave/solvers_maxCo_ddtCorrOff; fin=$wave/finalSolver_maxCo
i=0
for mesh in 95x20:0006:0007:0008 190x40:0018:0019:0020 380x80:0030:0031:0032; do
    IFS=: read -r m a b c <<< "$mesh"
    make1 "$off/$a" interFoam_$m
    make1 "$off/$b" interIsoFoam_$m
    make1 "$off/$c" interIsoGoam_$m
    make1 "$fin/$c" final_$m
done
