#!/bin/bash
# Campaign B (README.md): sampler on the desktop cores, driver and mpirun on cores 24-27.
cd "${0%/*}" || exit 1
source ~/OpenFOAM/OpenFOAM-v2512/etc/bashrc > /dev/null 2>&1
type blockMesh > /dev/null || { echo "OpenFOAM environment missing" >&2; exit 1; }
taskset -c 28-31 nice -n 19 python3 sampler.py log.sampler < /dev/null > /dev/null 2>&1 &
sampler=$!
taskset -c 24-27 python3 timing.py run "$@" > log.timing_B 2>&1 < /dev/null
kill $sampler
