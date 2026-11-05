"""Compares the finished final-solver cases with their Part 1 interIsoGoam counterparts (same case
number) using the metrics of analyseScan.py, and flags results that are worse than expected
from the earlier tests. Writes finalSolver_maxCo/check.txt. Usage: python3 finalSolver_checkMatrix.py

Flags: FAILED (solver stopped without 'End'); HEIGHT (at max Co <~ 0.2: |dH| more than 3
percentage points above max(|dH| of Part 1 interIsoGoam|, 2 %)); GROWTH (dH > +25 % at any Co);
UMAX (median max|U| more than 10 % above Part 1 interIsoGoam)."""
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
os.chdir(HERE)
from analyseScan import load, metrics  # noqa: E402

ref = load("solvers_fixedDt", "0036")
endTime = ref["t"][-1]
old = {}
for scan, f in [("finalSolver_fixedDt", "solvers_fixedDt.csv"), ("finalSolver_maxCo", "solvers_maxCo.csv")]:
    for r in csv.DictReader(open(f"data/{f}")):
        if r["application"] == "interIsoGoam":
            old[f"{scan}/{r['case']}"] = r
lines = [f"{'case':18s} {'mesh':>8s} {'dt/maxCo':>9s} {'Co':>5s} | {'dH':>6s} {'dH P1':>6s} | "
         f"{'Umed':>5s} {'P1':>5s} | {'steps':>7s} {'P1':>7s} | {'shape':>5s} {'P1':>5s} | flags"]
flags_all = []
for key in sorted(old, key=lambda k: (k.split('/')[0], int(old[k]['nx']), -float(old[k].get('deltaT') or old[k].get('maxCo')))):
    d = os.path.join(HERE, key)
    lg = os.path.join(d, "log.interIsoGoam")
    if not os.path.exists(lg):
        lg = os.path.join(d, "log.interIsoGoamExt")  # runs before the solver merge
    if not os.path.exists(lg):
        continue
    txt = open(lg, errors="ignore").read()
    o = old[key]
    par = o.get("deltaT") or o.get("maxCo")
    if "\nEnd\n" not in txt:
        t = [l for l in txt.splitlines() if l.startswith("Time = ")]
        alive = os.system(f"pgrep -f 'cd {d} ' > /dev/null") == 0
        state = "running" if alive else "STOPPED"
        if not alive:
            flags_all.append(f"{key}:FAILED")
        lines.append(f"{key:18s} {o['nx']+'x'+o['nz']:>8s} {par:>9s}  {state} at {t[-1] if t else '-'}")
        continue
    r = metrics(load(os.path.join(HERE, key.split('/')[0]), key.split('/')[1]), ref, endTime)
    co = float(r["avgMaxCo"]) if r.get("avgMaxCo") is not None else float("nan")
    dH, dH1 = r["dHfit"], float(o["dHfit"])
    U, U1 = r["Umed"], float(o["Umed"])
    fl = []
    if co <= 0.21 and abs(dH) > max(abs(dH1), 2.0) + 3.0: fl.append("HEIGHT")
    if dH > 25: fl.append("GROWTH")
    if U > 1.1*U1: fl.append("UMAX")
    flags_all += [f"{key}:{x}" for x in fl]
    lines.append(f"{key:18s} {o['nx']+'x'+o['nz']:>8s} {par:>9s} {co:5.2f} | {dH:+6.1f} {dH1:+6.1f} | "
                 f"{U:5.2f} {U1:5.2f} | {r['steps']:7d} {int(o['steps']):7d} | {r['shapeErr']:5.2f} "
                 f"{float(o['shapeErr']):5.2f} | {' '.join(fl)}")
open(os.path.join(HERE, "finalSolver_maxCo", "check.txt"), "w").write("\n".join(lines) + "\n")
print("\n".join(lines))
print("FLAGS: " + (", ".join(flags_all) if flags_all else "none"))
