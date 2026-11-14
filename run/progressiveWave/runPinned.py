# Copyright (C) 2026 Johan Roenby
# SPDX-License-Identifier: GPL-3.0-or-later
"""Runs the Part I progressive-wave studies on all cores of this machine, each run pinned.

The 760x160 runs use 8 processes on a pair of NUMA nodes (cores 0-7, 8-15, 16-23, 24-31),
longest first. The other runs are serial, one per core, on the last pair until they are
done; then that pair also takes 760x160 runs. Cases of solvers_maxCo, solvers_fixedDt and
nLayers_maxCo must have been created with parmScanner (without -exec); their output goes to
log.Allrun, so parmScanner -exec Allrun skips them afterwards. The ddtCorrOff studies are
built and run from these cases by ../common/runStoredDdtCorrOff.sh. A finished run leaves
.runPinned.done in its case and is skipped when the script is started again.
Steps per run, for the order, are taken from the archived runs of the same case number.
Needs the OpenFOAM v2512 environment. Usage: python3 runPinned.py [--dry-run]"""
import csv
import glob
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ARCHIVE = os.path.join(HERE, "..", "..", "..", "geometricGravity_archive", "rerun_20261111", "run",
                       "progressiveWave")
DICTATOR = os.path.expanduser("~/Projects/interFlow/dictator/dictator")
PAIRS = ["0-7", "8-15", "16-23", "24-31"]
SERIAL_PAIR = 3
MS_PER_STEP = {"95": 12, "190": 45, "380": 200, "760": 230}   # rough, for the order only


def old_steps(study, case):
    f = os.path.join(HERE, "data", f"{study}.csv")
    if os.path.exists(f):
        for r in csv.DictReader(open(f)):
            if r.get("case") == case and r.get("steps"):
                return int(r["steps"])
    return 20000


def parms(d):
    words = open(os.path.join(d, "parmScannerParms")).read().split()
    return {k.rstrip(":"): v for k, v in zip(words[1::2], words[2::2])}


def jobs():
    out = []
    for study in ["solvers_maxCo", "solvers_fixedDt", "nLayers_maxCo"]:
        for d in sorted(glob.glob(os.path.join(HERE, study, "[0-9][0-9][0-9][0-9]"))):
            p, c = parms(d), os.path.basename(d)
            nx = p.get("nx", "95")
            cost = old_steps(study, c) * MS_PER_STEP[nx]
            if p.get("application") == "streamFunFoam":
                out.append((d, "./Allrun > log.Allrun 2>&1", 1, 1))
                continue
            if nx == "760":
                cmd = ("dictator system/decomposeParDict numberOfSubdomains 8 > log.dictator && "
                       "dictator system/decomposeParDict method scotch >> log.dictator && "
                       "NUMA_LOCK_DIR=$PWD/.numaLock ./AllrunParallel > log.Allrun 2>&1 && "
                       "reconstructPar > log.reconstructPar 2>&1 && rm -rf processor*")
                out.append((d, cmd, 8, cost))
            else:
                out.append((d, "./Allrun > log.Allrun 2>&1", 1, cost))
            if study != "nLayers_maxCo":
                off = os.path.join(HERE, f"{study}_ddtCorrOff", c)
                np_ = 8 if nx == "760" else 1
                cmd = (f"mkdir -p {os.path.dirname(off)} && cd {HERE} && NUMA_LOCK_DIR={off}.numaLock "
                       f"../common/runStoredDdtCorrOff.sh {study}/{c} {study}_ddtCorrOff/{c} {np_} "
                       f"> {off}.log 2>&1 && ! grep -q FAILED {off}.log")
                out.append((off, cmd, np_, old_steps(f"{study}_ddtCorrOff", c) * MS_PER_STEP[nx]))
    return out


def done(d):
    return os.path.exists(os.path.join(d, ".runPinned.done"))


def main():
    dry = "--dry-run" in sys.argv
    todo = [j for j in jobs() if not done(j[0])]
    par = sorted([j for j in todo if j[2] == 8], key=lambda j: -j[3])
    ser = sorted([j for j in todo if j[2] == 1], key=lambda j: j[3])
    print(f"{len(par)} parallel runs, {len(ser)} serial runs; estimated "
          f"{sum(j[3] for j in par) / 3.6e6:.0f} h on 8 cores and {sum(j[3] for j in ser) / 3.6e6:.0f} h serial",
          flush=True)
    if dry:
        for j in par[:5] + ser[:5]:
            print(j)
        return
    first, last = map(int, PAIRS[SERIAL_PAIR].split("-"))
    free_pairs = [i for i in range(len(PAIRS)) if i != SERIAL_PAIR]
    free_cores = list(range(first, last + 1))
    running = {}   # Popen -> (job, resource)
    while par or ser or running:
        for proc, (job, res) in list(running.items()):
            if proc.poll() is None:
                continue
            del running[proc]
            status = "done" if proc.returncode == 0 else f"FAILED ({proc.returncode})"
            if proc.returncode == 0 and os.path.isdir(job[0]):
                open(os.path.join(job[0], ".runPinned.done"), "w").close()
            print(time.strftime("%F %T"), status, os.path.relpath(job[0], HERE), flush=True)
            if isinstance(res, int):
                free_cores.append(res)
            else:
                free_pairs.append(res[0])
        if not ser and not any(isinstance(r, int) for _, r in running.values()) and SERIAL_PAIR not in free_pairs \
                and SERIAL_PAIR not in [r[0] for _, r in running.values() if not isinstance(r, int)]:
            free_pairs.append(SERIAL_PAIR)
            free_cores = []
        while ser and free_cores:
            job, core = ser.pop(0), free_cores.pop(0)
            running[start(job, str(core))] = (job, core)
        while par and free_pairs:
            job, pair = par.pop(0), free_pairs.pop(0)
            running[start(job, PAIRS[pair])] = (job, (pair,))
        time.sleep(10)
    print(time.strftime("%F %T"), "all runs finished", flush=True)


def start(job, cpus):
    d, cmd = job[0], job[1]
    cwd = d if os.path.isdir(d) else HERE
    print(time.strftime("%F %T"), "start", os.path.relpath(d, HERE), "on cpus", cpus, flush=True)
    return subprocess.Popen(["taskset", "-c", cpus, "bash", "-c", f"source {DICTATOR} > /dev/null; {cmd}"],
                            cwd=cwd, stdin=subprocess.DEVNULL)


if __name__ == "__main__":
    main()
