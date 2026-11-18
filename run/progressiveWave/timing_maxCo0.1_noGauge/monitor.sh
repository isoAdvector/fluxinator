#!/bin/bash
# Copyright (C) 2026 Johan Roenby
# SPDX-License-Identifier: GPL-3.0-or-later
# Logs, every 10 s during the timing study: the load average, the clock frequency [MHz] of the
# first core of each NUMA node, and every process other than the timed runs that used more than
# 3 % of a core since the last sample (runPaired.sh writes this to log.monitor).
own="fluxinator|interIsoGoam|interFoam|interFoamDdtCorr|interIsoFoamDdtCorr|blockMesh|streamFunFoam|Allrun|runPaired|monitor|taskset|bash|sleep|cp|sed|awk|ps|cat|date"
while true; do
    freq=""
    for c in 0 4 8 12 16 20 24 28; do
        freq="$freq $(( $(cat /sys/devices/system/cpu/cpu$c/cpufreq/scaling_cur_freq) / 1000 ))"
    done
    echo "$(date +%T) load $(cut -d' ' -f1-3 /proc/loadavg) MHz$freq"
    # The second frame of top: CPU use over the 2 s between the frames
    top -b -n 2 -d 2 -w 200 | awk -v own="^($own)$" '
        /^top -/ {frame++} frame == 2 && seen && $9 > 3 && $12 !~ own {print "    other:", $12, $9 "%"}
        frame == 2 && /PID/ {seen = 1}'
    sleep 8
done
