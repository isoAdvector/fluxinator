# Copyright (C) 2026 Johan Roenby
# SPDX-License-Identifier: GPL-3.0-or-later
"""Runs the fluxinator matrix (finalSolver_fixedDt, finalSolver_maxCo) on at most MAXCORES cores, cheapest
cases first. Cost estimate: clock time of the Part 1 interIsoGoam run with the same case number
(data/solvers_*.csv). 760x160 cases run with AllrunParallel on 4 cores (one NUMA node) and
are reconstructed afterwards; all others serially with Allrun. A case counts as finished when its
solver log ends with 'End'.
Progress: log.matrix. Usage: nohup python3 finalSolver_runMatrix.py > log.runMatrix 2>&1 &"""
import csv
import os
import subprocess
import time

MAXCORES = 16
HERE = os.path.dirname(os.path.abspath(__file__))
ENV = ""


def log(msg):
    with open(os.path.join(HERE, "log.matrix"), "a") as f:
        f.write(time.strftime("%Y-%m-%d %H:%M:%S ") + msg + "\n")


cost = {}
for scan, f in [("finalSolver_fixedDt", "solvers_fixedDt.csv"), ("finalSolver_maxCo", "solvers_maxCo.csv")]:
    for r in csv.DictReader(open(os.path.join(HERE, "data", f))):
        if r["application"] == "interIsoGoam":
            cost[f"{scan}/{r['case']}"] = (float(r["clockTime"]), int(r["nx"]))


def finished(case):
    for app in ("fluxinator", "interIsoGoam", "interIsoGoamExt"):  # the last two: runs before fluxinator
        lg = os.path.join(HERE, case, "log." + app)
        if os.path.exists(lg):
            return "\nEnd\n" in open(lg, errors="ignore").read()
    return False


cost = {c: v for c, v in cost.items() if os.path.isdir(os.path.join(HERE, c))}
jobs = sorted([c for c in cost if not finished(c)], key=lambda c: cost[c][0])
log(f"start: {len(jobs)} cases to run, {len(cost) - len(jobs)} already finished")
running = {}
while jobs or running:
    for c, (p, n) in list(running.items()):
        if p.poll() is not None:
            del running[c]
            log(f"done  {c} (exit {p.returncode}, {'End' if finished(c) else 'NO End'})")
    used = sum(n for _, n in running.values())
    while jobs:
        c = jobs[0]; n = 4 if cost[c][1] == 760 else 1
        if used + n > MAXCORES:
            break
        jobs.pop(0)
        d = os.path.join(HERE, c)
        if n == 1:
            cmd = ENV + f"cd {d} && ./Allrun > log.Allrun 2>&1"
        else:
            cmd = ENV + (f"cd {d} && ./AllrunParallel > log.Allrun 2>&1; "
                         f"reconstructPar > log.reconstructPar 2>&1 && rm -rf processor*")
        running[c] = (subprocess.Popen(["bash", "-c", cmd]), n)
        used += n
        log(f"start {c} on {n} core(s), est. {cost[c][0]/3600*(1 if n == 1 else 0.3):.2f} h")
    time.sleep(30)
log("all cases finished")
