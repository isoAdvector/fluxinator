#!/bin/bash
# Copyright (C) 2026 Johan Roenby
# SPDX-License-Identifier: GPL-3.0-or-later
# Timing study at max Co 0.1 without the wave gauge (surfacePosition1 enabled false in
# controlDict.run; made when the gauge was interfaceHeight, which rebuilds a mesh search
# tree every step; the gauge is now surfacePosition, whose cost is negligible). 36 runs to
# t = 6 s, largest meshes first, each pinned with taskset to one of the cores in CORES (by
# default the first core of each 4-core NUMA node of the 32-core machine, so that each run
# has its own memory controller and L3 cache; memory is then allocated on that node). The
# runs are dealt to the cores in turn and each core runs its share one after the other.
# Nothing else should run on these NUMA nodes. Needs the OpenFOAM environment of the runs
# (v2512) and dictator. Usage: [CORES="0 4 ..."] ./runAll.sh
cd "${0%/*}" || exit
[ -d final_95x20_r1 ] || ./makeCases.sh || exit 1
read -r -a cores <<< "${CORES:-0 4 8 12 16 20 24 28}"
cases=(*_380x80_r? *_190x40_r? *_95x20_r?)
for i in "${!cores[@]}"; do
    (
        for ((j = i; j < ${#cases[@]}; j += ${#cores[@]})); do
            (cd "${cases[j]}" && taskset -c "${cores[i]}" ./Allrun > /dev/null 2>&1)
        done
    ) &
done
wait
