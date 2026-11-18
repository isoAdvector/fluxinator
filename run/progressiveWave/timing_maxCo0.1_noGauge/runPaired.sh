#!/bin/bash
# Copyright (C) 2026 Johan Roenby
# SPDX-License-Identifier: GPL-3.0-or-later
# Timing study at max Co 0.1 in paired rounds: in each round the five solvers of one mesh run at
# the same time, each pinned with taskset to the first core of its own 4-core NUMA node (memory
# is then allocated on that node), so that a slower or faster state of the machine affects all
# five alike and cancels in their ratios. Three rounds per mesh, largest mesh first; the nodes
# are rotated from round to round. monitor.sh logs anything else that uses the machine
# (log.monitor). Needs the OpenFOAM environment of the runs (v2512) and dictator.
# Usage: ./runPaired.sh   (then python3 analyse.py)
cd "${0%/*}" || exit
[ -d final_95x20_r1 ] || ./makeCases.sh || exit 1
cores=(0 4 8 12 16 20 24 28)
solvers=(final interFoamDefault interFoam interIsoFoam interIsoGoam)
./monitor.sh > log.monitor 2>&1 &
monitor=$!
g=0
for mesh in 380x80 190x40 95x20; do
    for r in 1 2 3; do
        echo "$(date +%T) round $g: $mesh repeat $r" >> log.rounds
        for k in "${!solvers[@]}"; do
            c=${cores[$(( (k + 5*g) % ${#cores[@]} ))]}
            d=${solvers[k]}_${mesh}_r$r
            echo "    $d on core $c" >> log.rounds
            (cd "$d" && taskset -c "$c" ./Allrun > /dev/null 2>&1) &
        done
        wait $(jobs -p | grep -v "^$monitor$")
        g=$((g + 1))
    done
done
echo "$(date +%T) done" >> log.rounds
kill $monitor
