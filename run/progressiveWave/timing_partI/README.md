<!--
  Copyright (C) 2026 Johan Roenby
  SPDX-License-Identifier: GPL-3.0-or-later
-->
# Timing of the Part I progressive-wave runs (Fig. waveCost b)

`timing.py make | run | analyse` (OpenFOAM v2512, dictator). The design is described in the
docstring of `timing.py`: copies of the `../solvers_maxCo` set-ups with a shorter `endTime`,
time per step over t > 0.25 s, compared solvers run at the same time, two repeats with the
solvers on other cores. `analyse` writes `../data/timing_partI.csv`; the figure multiplies the
time per step by the number of steps of the full run (`papers/part1/fig_waveCost`).

## Conditions of the run on 2026-11-13, 10:14 onwards

- 2x AMD EPYC 7351 (32 cores, 8 NUMA nodes of 4 cores), governor `performance`, boost off
  (2.4 GHz), `vm.zone_reclaim_mode` 1. The page cache was not dropped; before the start every
  node had 4-14 GB free (43 GB cache in total).
- The Part I rerun (`../runPinned.py`) was paused with SIGSTOP: its 41 processes held about
  9 GB of memory but used no CPU. Otherwise only the desktop session was running.
- Binding: every solver rank on its own core with a rankfile and
  `hwloc_base_mem_alloc_policy local_only`, also the serial runs (`mpirun -np 1`); pre-processing
  pinned with taskset. `mpirun --cpu-set ... --bind-to core` was not used: it binds every rank
  to the whole set. `log.bindings` in each case holds the report of `--report-bindings`.
- Checked during round 0 (/proc/<pid>/numa_maps): anonymous memory of every solver 100 % on
  its own node; only file-backed pages (shared libraries, ~140 MB) were mostly remote.
- Steadiness: `msFirstHalf` and `msSecondHalf` in the csv give the time per step in the two
  halves of the measured window; they agree within 5 % in every run.
- Disturbance (campaign A, `log.timing`, `log.timing_r3`): the 760x160 rounds used all 32
  cores, so every other process had to share a core with a solver rank, and an 8-rank job
  waits for its slowest rank. In r1 the rounds at maxCo 0.5 and 0.2 (10:33-10:41) gave
  interIsoGoam on cores 16-23 315 and 290 ms/step, against 252 and 245 in r2; in r3
  (11:15-11:45) all rounds were unsteady (halves differing by 10-40 %), during a browser
  sign-in and an automatic update of a VS Code extension (11:41, journal). The serial rounds
  left 8 cores free and were steady. The 760x160 runs of campaign A are kept in
  `disturbed760/` and not used.

## Campaign B, 2026-11-13 from 12:05 (`runB.sh`, `log.timing_B`, `log.sampler`)

- All running processes of the user (130: desktop, ThinLinc, editor, browser, Claude
  clients; not the paused rerun) confined to cores 28-31 with `taskset -a -p` (previous
  affinities in `affinity_before.txt`; restored afterwards); the timing driver and mpirun on
  cores 24-27; processes of other users (root) not confined.
- 760x160: three 8-rank jobs per round on cores 0-23, repeats r4-r6 with the solvers rotated
  over the pairs of nodes. Second layer: interIsoGoam nLayers 2 and 1 on 380x80 (r3, r4).
- `sampler.py` (cores 28-31, nice 19) logs every 5 s each other thread above 2 % of a core
  with its core, and the sum of such use on cores 0-23 (`SUM` lines).
- A first launch at 12:05:05 failed within seconds (OpenFOAM environment missing; logs
  `*_failedLaunch`).
