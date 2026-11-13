# Copyright (C) 2026 Johan Roenby
# SPDX-License-Identifier: GPL-3.0-or-later
"""Logs CPU use outside the timing runs: every INTERVAL s, every thread (all users, also
kernel threads) that is not a solver or mpirun and used more than THRESHOLD % of a core,
with the core it last ran on; and per sample the summed use of such threads on the timing
cores. Reads only /proc. Run it on the cores kept free for the desktop, e.g.
taskset -c 28-31 nice -n 19 python3 sampler.py log.sampler"""
import glob
import os
import sys
import time

INTERVAL = 5
THRESHOLD = 2.0
TIMING_CORES = set(range(0, 24))
SKIP = {"interFoam", "interIsoFoam", "interIsoGoam", "mpirun", "orted", "blockMesh",
        "streamFunFoam", "decomposePar"}
TICK = os.sysconf("SC_CLK_TCK")


def sample():
    out = {}
    for f in glob.glob("/proc/[0-9]*/task/[0-9]*/stat"):
        try:
            s = open(f).read()
        except OSError:
            continue
        comm = s[s.index("(") + 1:s.rindex(")")]
        w = s[s.rindex(")") + 2:].split()
        # after the comm: state is field 3, utime 14, stime 15, processor 39 (1-based)
        out[f.split("/")[4]] = (f.split("/")[2], comm, int(w[11]) + int(w[12]), int(w[36]))
    return out


def main():
    log = open(sys.argv[1], "a", buffering=1)
    log.write(f"# {time.strftime('%F %T')} start; threads > {THRESHOLD} % over {INTERVAL} s: "
              "time pid tid comm %cpu core; SUM: time %cpu of such threads on cores 0-23\n")
    old, t_old = sample(), time.time()
    while True:
        time.sleep(INTERVAL)
        new, t = sample(), time.time()
        busy = 0.0
        stamp = time.strftime("%T")
        for tid, (pid, comm, ticks, core) in new.items():
            if comm in SKIP or tid not in old:
                continue
            pct = 100 * (ticks - old[tid][2]) / TICK / (t - t_old)
            if pct > THRESHOLD:
                log.write(f"{stamp} {pid} {tid} {comm.replace(' ', '_')} {pct:.0f} {core}\n")
                if core in TIMING_CORES:
                    busy += pct
        log.write(f"SUM {stamp} {busy:.0f}\n")
        old, t_old = new, t


if __name__ == "__main__":
    main()
