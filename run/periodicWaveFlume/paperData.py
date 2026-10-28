#!/usr/bin/env python3
#
# Writes the numbers behind the wave results in the paper to paperData/
# (tracked in git, unlike the scans and the metrics files):
#   runs_fixedDt.csv, runs_maxCo.csv       copies of the metrics written by
#                                          analyseScan.py (see there for the
#                                          definitions)
#   statistics.txt                         derived numbers quoted in the text
#
# Usage (from this directory, after Allrun has written fixedDtMetrics.csv and
# maxCoMetrics.csv):
#   python3 paperData.py

import csv
import math
import os
import re
import shutil
import statistics

import numpy as np

OUT = "paperData"
FIXED_CSV = "fixedDtMetrics.csv"
MAXCO_CSV = "maxCoMetrics.csv"
FIXED_DIR = "gravityModelScan_fixedDt"
MAXCO_DIR = "gravityModelScan"
REF = "0036"
U_CREST = 0.44453480326804
SOLVERS = ["interFoam", "interIsoFoam", "interIsoGoam"]
NZ = [20, 40, 80, 160]
NX = {20: 95, 40: 190, 80: 380, 160: 760}


def load(path, key):
    rows = []
    with open(path) as f:
        for r in csv.DictReader(f):
            d = {}
            for k, v in r.items():
                if k in ("case", "application"):
                    d[k] = v
                elif k == "diverged":
                    d[k] = (v == "True")
                else:
                    try:
                        d[k] = float(v)
                    except ValueError:
                        d[k] = None
            d["nz"] = int(d["nz"])
            d["level"] = d[key]
            if d["application"] in SOLVERS:
                rows.append(d)
    return rows


def run_of(rows, app, nz, i):
    """i-th run of a solver on a mesh, ordered from the largest dt or maxCo."""
    return sorted((r for r in rows if r["application"] == app and r["nz"] == nz),
                  key=lambda r: -r["level"])[i]


def fmt(v, spec, missing="--"):
    return missing if v is None else format(v, spec)


def interp_fixed(fix, r):
    pts = sorted((math.log(q["deltaT"]), q["dHfit"]) for q in fix
                 if q["application"] == r["application"] and q["nz"] == r["nz"]
                 and q["dHfit"] is not None)
    x = math.log(30.0/r["steps"])
    if x < pts[0][0] or x > pts[-1][0]:
        return None
    for (x0, y0), (x1, y1) in zip(pts, pts[1:]):
        if x0 <= x <= x1:
            return y0 + (y1 - y0)*(x - x0)/(x1 - x0)
    return None


def reference_wave():
    txt = open(os.path.join(FIXED_DIR, REF, "constant", "out")).read()
    val = {k: float(re.search(rf"^{k}\s+(\S+);", txt, re.M).group(1))
           for k in ("H", "d", "T", "lambda", "k", "h_crest", "h_trough", "u_crest")}
    t, h = np.loadtxt(os.path.join(FIXED_DIR, REF, "postProcessing", "interfaceHeight1", "0",
                                   "height.dat"), comments="#", usecols=(0, 2), unpack=True)
    return val, h.max(), h.min()


def max_u_series(case):
    return np.loadtxt(os.path.join(FIXED_DIR, case, "postProcessing", "fieldMinMax1", "0",
                                   "fieldMinMax.dat"), comments="#", usecols=(0, 2), unpack=True)


def main():
    os.makedirs(OUT, exist_ok=True)
    shutil.copy(FIXED_CSV, os.path.join(OUT, "runs_fixedDt.csv"))
    shutil.copy(MAXCO_CSV, os.path.join(OUT, "runs_maxCo.csv"))
    F = load(FIXED_CSV, "deltaT")
    M = load(MAXCO_CSV, "maxCo")

    lines = []
    p = lines.append
    campaigns = (("fixed dt", F), ("adaptive", M))

    val, crest, trough = reference_wave()
    p("Exact wave (streamFunctionWave output, " + f"{FIXED_DIR}/{REF}/constant/out)")
    p(f"  H = {val['H']} m, d = {val['d']} m, T = {val['T']} s, L = {val['lambda']:.4f} m, "
      f"kd = {val['k']*val['d']:.3f}, H/d = {val['H']/val['d']:.4f}, H/L = {val['H']/val['lambda']:.4f}")
    p(f"  crest {1e3*val['h_crest']:.2f} mm, trough {1e3*val['h_trough']:.2f} mm, "
      f"u_crest {val['u_crest']:.4f} m/s")
    p(f"  sampled at the gauge in {REF}: crest {1e3*crest:.2f} mm, trough {1e3*trough:.2f} mm, "
      f"H {1e3*(crest - trough):.3f} mm")

    p("\nAverage max Courant number at the largest and smallest fixed time step")
    for i in (0, 3):
        cos = [run_of(F, a, n, i)["avgMaxCo"] for a in SOLVERS for n in NZ]
        cos = [c for c in cos if c is not None]
        p(f"  dt0/{2**i}: {min(cos):.3f} ... {max(cos):.3f}")

    p("\nWave height error: local linear trend at the last wave minus height of the last wave")
    for name, R in campaigns:
        d = [abs(r["dHfit"] - r["dHend"]) for r in R if r["dHfit"] is not None]
        p(f"  {name}: median {statistics.median(d):.2f}, max {max(d):.2f} percentage points ({len(d)} runs)")

    p("\n380x80 mesh, fixed dt: trend and last-wave height error, last crest and trough, phase, Courant number")
    for i in (0, 2):
        for a in SOLVERS:
            r = run_of(F, a, 80, i)
            p(f"  dt = {1e3*r['deltaT']:g} ms {a:12s} dHfit {r['dHfit']:+.1f} %  dHend {r['dHend']:+.1f} %  "
              f"crest {r['crestEnd']:+.1f} %  trough {r['troughEnd']:+.1f} %  phase {r['phase']:+.3f} s  "
              f"Co {r['avgMaxCo']:.3f}")

    p("\nRatio of the wave height errors of interIsoFoam and interIsoGoam at the smallest dt / maxCo")
    for name, R in campaigns:
        p(f"  {name}: " + ", ".join(
            f"{run_of(R, 'interIsoFoam', n, 3)['dHfit']/run_of(R, 'interIsoGoam', n, 3)['dHfit']:.2f}"
            for n in NZ))

    p("\nMean celerity error at the smallest dt / maxCo, all solvers and meshes")
    for name, R in campaigns:
        dc = [run_of(R, a, n, 3)["dcMean"] for a in SOLVERS for n in NZ]
        dc = [v for v in dc if v is not None]
        p(f"  {name}: {min(dc):+.3f} ... {max(dc):+.3f} % (max |dc/c| = {max(map(abs, dc)):.3f} %)")
    p("Mean celerity error at the largest dt / maxCo, and the spread across solvers")
    for name, R in campaigns:
        for n in NZ:
            v = [run_of(R, a, n, 0)["dcMean"] for a in SOLVERS]
            got = [x for x in v if x is not None]
            spread = f"{max(got)-min(got):.3f}" if len(got) == 3 else "--"
            p(f"  {name} {NX[n]}x{n}: " + ", ".join(
                f"{a} {fmt(x, '+.3f')}" for a, x in zip(SOLVERS, v)) + f"  spread {spread}")
    p("Crest arrival-time (phase) error at the smallest dt / maxCo, all solvers and meshes")
    for name, R in campaigns:
        ph = [run_of(R, a, n, 3)["phase"] for a in SOLVERS for n in NZ]
        p(f"  {name}: {min(ph):+.3f} ... {max(ph):+.3f} s")
    p("Crest arrival-time (phase) error at the largest dt / maxCo")
    for name, R in campaigns:
        for n in NZ:
            p(f"  {name} {NX[n]}x{n}: " + ", ".join(
                f"{a} {fmt(run_of(R, a, n, 0)['phase'], '+.3f')} (Co {fmt(run_of(R, a, n, 0)['avgMaxCo'], '.2f')})"
                for a in SOLVERS))

    p("\nRuns whose crests carry more than one peak (mean peaks per crest, see analyseScan.py)")
    for name, R in campaigns:
        for r in R:
            if r["crestPeaks"] is not None and r["crestPeaks"] != 1:
                p(f"  {name} {r['case']} {r['application']} {NX[r['nz']]}x{r['nz']} level {r['level']:g}: "
                  f"{r['crestPeaks']:.2f}")

    p("\nMedian max|U|/u_crest: runs where interIsoGoam is not the lowest")
    bad = [(name, n, i) for name, R in campaigns for n in NZ for i in range(4)
           if any(run_of(R, a, n, i)["Umed"] is not None
                  and run_of(R, a, n, i)["Umed"] <= run_of(R, "interIsoGoam", n, i)["Umed"]
                  for a in SOLVERS[:2])]
    p(f"  {bad if bad else 'none'}")

    p("\nmax|U|/u_crest over time, 760x160 mesh, dt0/8 (fixed dt)")
    for case in ("0046", "0047", "0048"):
        t, u = max_u_series(case)
        r = u/U_CREST
        late = t > 2
        p(f"  {case}: max {r.max():.3f}, t > 2 s: {r[late].min():.3f} ... {r[late].max():.3f}, "
          f"median 2-7 s {np.median(r[late & (t < 7)]):.3f}, median 7-30 s {np.median(r[t >= 7]):.3f}, "
          f"10-90 % percentiles t > 10 s {np.percentile(r[t > 10], 10):.2f} ... "
          f"{np.percentile(r[t > 10], 90):.2f}, final {r[-1]:.3f}")

    p("\nDiverged runs")
    for name, R in campaigns:
        for r in R:
            if r["diverged"]:
                p(f"  {name} {r['case']} {r['application']} {NX[r['nz']]}x{r['nz']} level {r['level']:g} "
                  f"stopped at t = {r['end']:g} s")

    p("\nWall time per step [ms], fixed dt, all solvers and time steps")
    for n in NZ:
        ms = [run_of(F, a, n, i)["msPerStep"] for a in SOLVERS for i in range(4)]
        p(f"  {NX[n]}x{n}: {min(ms):.1f} ... {max(ms):.1f}")

    p("\nChange in the number of time steps of interIsoGoam relative to interFoam / interIsoFoam [%]")
    for n in NZ:
        p(f"  {NX[n]}x{n}: " + ", ".join(
            f"maxCo {run_of(M, 'interIsoGoam', n, i)['maxCo']:g}: "
            f"{100*(run_of(M, 'interIsoGoam', n, i)['steps']/run_of(M, 'interFoam', n, i)['steps'] - 1):+.0f} / "
            f"{100*(run_of(M, 'interIsoGoam', n, i)['steps']/run_of(M, 'interIsoFoam', n, i)['steps'] - 1):+.0f}"
            for i in range(4)))
    p("  interFoam relative to interIsoFoam, 760x160, maxCo 0.1 and 0.05: " + ", ".join(
        f"{100*(run_of(M, 'interFoam', 160, i)['steps']/run_of(M, 'interIsoFoam', 160, i)['steps'] - 1):+.0f}"
        for i in (2, 3)))

    p("\nRuns with |dH| <= 2 % and shape error <= 1 %")
    for name, R in campaigns:
        for a in SOLVERS:
            ok = [r for r in R if r["application"] == a and r["dHfit"] is not None
                  and abs(r["dHfit"]) <= 2 and r["shapeErr"] <= 1]
            p(f"  {name} {a}: {len(ok)} runs: " + ", ".join(
                f"{NX[r['nz']]}x{r['nz']} level {r['level']:g} (Co {r['avgMaxCo']:.2f})" for r in ok))

    p("\nAdaptive runs against the fixed-dt trend interpolated to their mean time step")
    diffs, exceptions = [], []
    for r in M:
        if r["dHfit"] is None:
            continue
        yi = interp_fixed(F, r)
        if yi is None:
            continue
        diffs.append(abs(r["dHfit"] - yi))
        if abs(r["dHfit"] - yi) > 2:
            exceptions.append(f"{r['application']} {NX[r['nz']]}x{r['nz']} maxCo {r['maxCo']:g}: "
                              f"{r['dHfit']:+.1f} % against {yi:+.1f} %")
    p(f"  {len(diffs)} of {len(M)} runs comparable, median difference {statistics.median(diffs):.2f} "
      f"percentage points, {sum(d <= 2 for d in diffs)} within 2 points")
    for e in exceptions:
        p("  exception: " + e)

    with open(os.path.join(OUT, "statistics.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")
    print("Wrote", ", ".join(sorted(os.listdir(OUT))))


if __name__ == "__main__":
    main()
