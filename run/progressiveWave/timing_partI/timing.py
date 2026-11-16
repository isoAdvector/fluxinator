# Copyright (C) 2026 Johan Roenby
# SPDX-License-Identifier: GPL-3.0-or-later
"""Timing of the Part I progressive-wave runs under a Courant limit (Fig. waveCost b).

The wall-clock time of a full run is its number of time steps (the full run in
../solvers_maxCo) times the wall-clock time per step, measured here in short runs of the
same set-ups on an otherwise idle machine. The time per step is steady from about 0.25 s
of simulated time on, so it is measured over t > 0.25 s.

Runs that are compared run at the same time, each solver bound to its own cores with local
memory (mpirun with a rankfile, one core per rank, and hwloc_base_mem_alloc_policy local_only;
mpirun --cpu-set binds every rank to the whole set), also the serial ones:
- 95x20, 190x40, 380x80: the three solvers of one mesh and Courant limit on three cores of
  one NUMA node (the fourth core idle), eight such groups at a time;
- 760x160: interFoam, interIsoFoam and interIsoGoam, 8 ranks each on two NUMA nodes (cores
  0-23), one Courant limit at a time, repeats r4, r5, r6 with the solvers rotated over the
  three pairs of nodes; all other processes of the user (desktop, editor, browser) confined
  to cores 28-31 with taskset, the driver and mpirun on cores 24-27. The first 760x160 runs
  (r1-r3, all 32 cores, desktop not confined) were disturbed by other processes and are kept
  in disturbed760/ (see README.md);
- the cost of the second layer: interIsoGoam with nLayers 2 and 1 on 380x80, repeats r3, r4,
  two such pairs per NUMA node 0-3.
When the first run of a group ends, the others of the group are stopped, so that all
measured steps ran under the same load. The coarser meshes are repeated once (r2) with the
solvers on other cores. Each case is a copy of the solvers_maxCo set-up with a shorter
endTime. sampler.py logs the CPU use of all other threads during the runs.

Needs the OpenFOAM v2512 environment and dictator.
Usage: python3 timing.py make | run [round ...] | analyse   (make keeps existing cases;
run without round numbers runs all rounds)"""
import csv
import glob
import os
import shutil
import signal
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "..", "solvers_maxCo")
DATA = os.path.join(HERE, "..", "data")
DICTATOR = os.path.expanduser("~/Projects/interFlow/dictator/dictator")
END = {"95": 30, "190": 4, "380": 2, "760": {"0.5": 2, "0.2": 2, "0.1": 1, "0.05": 1}}
T_START = 0.25
SOLVERS = ["interFoam", "interIsoFoam", "interIsoGoam"]
COS = ["0.5", "0.2", "0.1", "0.05"]


def parms(d):
    words = open(os.path.join(d, "parmScannerParms")).read().split()
    return {k.rstrip(":"): v for k, v in zip(words[1::2], words[2::2])}


def source_cases():
    """(nx, maxCo, application) -> solvers_maxCo case"""
    out = {}
    for d in sorted(glob.glob(os.path.join(SRC, "[0-9][0-9][0-9][0-9]"))):
        p = parms(d)
        if p["application"] in SOLVERS:
            out[(p["nx"], p["maxCo"], p["application"])] = os.path.basename(d)
    return out


def name(nx, co, app, nL=2):
    return f"{app}{'_nL1' if nL == 1 else ''}_{nx}_Co{co}"


def make():
    src = source_cases()
    jobs = [(nx, co, app, 2, (4, 5, 6) if nx == "760" else (1, 2)) for (nx, co, app) in src]
    jobs += [("380", co, "interIsoGoam", nL, (3, 4)) for co in COS for nL in (2, 1)]
    n = 0
    for nx, co, app, nL, repeats in jobs:
        for r in repeats:
            d = os.path.join(HERE, f"{name(nx, co, app, nL)}_r{r}")
            if os.path.exists(d):
                continue
            n += 1
            s = os.path.join(SRC, src[(nx, co, app)])
            os.makedirs(d)
            for f in ["0.orig", "system"]:
                shutil.copytree(os.path.join(s, f), os.path.join(d, f))
            shutil.copytree(os.path.join(s, "constant"), os.path.join(d, "constant"),
                            ignore=shutil.ignore_patterns("polyMesh"))
            end = END[nx][co] if nx == "760" else END[nx]
            cmds = [f"dictator system/controlDict.run endTime {end}"]
            if nx == "760":
                cmds += ["dictator system/decomposeParDict numberOfSubdomains 8",
                         "dictator system/decomposeParDict method scotch"]
            if nL == 1:
                cmds += ["dictator system/fvSchemes interfaceForces.gravity.nLayers 1"]
            subprocess.run(["bash", "-c", f"source {DICTATOR} > /dev/null; set -e; " + "; ".join(cmds)],
                           cwd=d, check=True, stdout=subprocess.DEVNULL)
            open(os.path.join(d, "timing.parms"), "w").write(
                f"source {src[(nx, co, app)]}\napplication {app}\nnx {nx}\nmaxCo {co}\nnLayers {nL}\n"
                f"repeat {r}\n")
    print(f"{n} new cases")


def groups():
    """Rounds of groups; a group is a list of (case, cpus)"""
    rounds = []
    for r in (1, 2):
        # 380x80 and 190x40: one mesh and Courant limit per NUMA node; r2 on the other nodes,
        # with the solvers rotated over the three cores
        order = SOLVERS if r == 1 else SOLVERS[1:] + SOLVERS[:1]
        triples = [(nx, co) for nx in ("380", "190") for co in COS]
        if r == 2:
            triples = triples[4:] + triples[:4]
        rounds.append([[(f"{name(nx, co, app)}_r{r}", str(4 * node + k)) for k, app in enumerate(order)]
                       for node, (nx, co) in enumerate(triples)])
    rounds.append([[(f"{name('95', co, app)}_r{r}", str(4 * node + k))
                    for k, app in enumerate(SOLVERS if r == 1 else SOLVERS[1:] + SOLVERS[:1])]
                   for r in (1, 2) for node, co in enumerate(COS, start=0 if r == 1 else 4)])
    pairs = ["0-7", "8-15", "16-23"]
    for r in (4, 5, 6):                                   # rounds 3-6, 7-10, 11-14
        for co in COS:
            rounds.append([[(f"{name('760', co, a)}_r{r}", pairs[(k + r - 4) % 3])
                            for k, a in enumerate(SOLVERS)]])
    # Second layer, 380x80 (round 15): nLayers 2 and 1 side by side, r4 with the cores swapped
    rounds.append([[(f"{name('380', co, 'interIsoGoam', nL)}_r{r}", str(2 * g + k))
                    for k, nL in enumerate((2, 1) if r == 3 else (1, 2))]
                   for g, (r, co) in enumerate((r, co) for r in (3, 4) for co in COS)])
    return rounds


def expand(cpus):
    a, _, b = cpus.partition("-")
    return ",".join(str(c) for c in range(int(a), int(b or a) + 1))


def start(case, cpus):
    d = os.path.join(HERE, case)
    p = dict(l.split() for l in open(os.path.join(d, "timing.parms")))
    app, n = p["application"], len(expand(cpus).split(","))
    pre = (f"rm -rf 0 processor* && cp -r 0.orig 0 && touch periodicWaveTank.foam && "
           f"cp system/controlDict.init system/controlDict && "
           f"taskset -c {cpus} blockMesh > log.blockMesh 2>&1 && "
           f"taskset -c {cpus} streamFunFoam > log.streamFunFoam 2>&1 && "
           f"cp 0.orig/p_rgh 0/p_rgh && cp system/controlDict.run system/controlDict")
    if n > 1:
        pre += f" && taskset -c {cpus} decomposePar > log.decomposePar 2>&1"
    open(os.path.join(d, "rankfile"), "w").write(
        "".join(f"rank {i}=localhost slot={c}\n" for i, c in enumerate(expand(cpus).split(","))))
    run = (f"mpirun --rankfile rankfile --mca rmaps_rank_file_physical true --report-bindings "
           f"--mca hwloc_base_mem_alloc_policy local_only -np {n} {app}{' -parallel' if n > 1 else ''} "
           f"< /dev/null > log.{app} 2> log.bindings")
    open(os.path.join(d, "timing.cpus"), "w").write(cpus + "\n")
    print(time.strftime("%F %T"), "start", case, "on cpus", cpus, flush=True)
    return subprocess.Popen(["bash", "-c", f"{pre} && {run}"], cwd=d, stdin=subprocess.DEVNULL,
                            start_new_session=True)


def run(sel=None):
    for i, rnd in enumerate(groups()):
        if sel and i not in sel:
            continue
        print(time.strftime("%F %T"), f"round {i}", flush=True)
        live = [{c: start(c, cpus) for c, cpus in g} for g in rnd]
        while any(live):
            time.sleep(5)
            for g in live:
                ended = [c for c, p in g.items() if p.poll() is not None]
                if not ended:
                    continue
                for c in ended:
                    print(time.strftime("%F %T"), "ended" if g[c].returncode == 0 else
                          f"FAILED ({g[c].returncode})", c, flush=True)
                    del g[c]
                for c, p in g.items():   # stop the rest of the group
                    os.killpg(p.pid, signal.SIGTERM)
                    p.wait()
                    print(time.strftime("%F %T"), "stopped", c, flush=True)
                g.clear()
    print(time.strftime("%F %T"), "all rounds finished", flush=True)


def ms_per_step(log):
    t, T, E = None, [], []
    for line in open(log):
        if line.startswith("Time = "):
            t = float(line.split()[2])
        elif line.startswith("ExecutionTime") and t is not None:
            w = line.split()
            T.append(t)
            E.append(float(w[2]))
    idx = [i for i, x in enumerate(T) if x > T_START]
    n = len(idx) - 1
    if n < 2:
        return float("nan"), n, T[-1] if T else 0, float("nan"), float("nan")
    a, m, b = idx[0], idx[n // 2], idx[-1]
    return (1000 * (E[b] - E[a]) / n, n, T[-1], 1000 * (E[m] - E[a]) / (m - a),
            1000 * (E[b] - E[m]) / (b - m))


def analyse():
    rows = []
    for f in sorted(glob.glob(os.path.join(HERE, "*_r[1-6]", "timing.parms"))):
        d = os.path.dirname(f)
        p = dict(l.split() for l in open(f))
        ms, n, tEnd, ms1, ms2 = ms_per_step(os.path.join(d, f"log.{p['application']}"))
        rows.append({**p, "case": os.path.basename(d), "cpus": open(os.path.join(d, "timing.cpus")).read().strip(),
                     "msPerStep": f"{ms:.3f}", "stepsMeasured": n, "tEnd": tEnd,
                     "msFirstHalf": f"{ms1:.3f}", "msSecondHalf": f"{ms2:.3f}"})
    out = os.path.join(DATA, "timing_partI.csv")
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"{len(rows)} runs -> {os.path.relpath(out, HERE)}")


if __name__ == "__main__":
    if sys.argv[1] == "run":
        run([int(a) for a in sys.argv[2:]])
    else:
        {"make": make, "analyse": analyse}[sys.argv[1]]()
