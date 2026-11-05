"""Pressure-solver study: accuracy metrics (as papers/part2/finalSolver.py) and cost
of each variant. Writes summary.txt. Run with python3 analyse.py."""
import glob, os, re, sys
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "../../../papers/part2"))
import finalSolver as F  # noqa: E402

REF = F.load(os.path.join(F.WAVE, "solvers_fixedDt"), "0036")


def cost(d):
    log = [l for l in glob.glob(d + "/log.*") if not re.search("blockMesh|streamFun", l)][0]
    steps, it, t, ex = 0, 0, 0.0, 0.0
    for line in open(log):
        if line.startswith("Time = "):
            steps += 1
            t = float(line.split()[2])
        elif "Solving for p_rgh" in line:
            it += int(line.rsplit(" ", 1)[1])
        elif line.startswith("ExecutionTime"):
            ex = float(line.split()[2])
    return dict(steps=steps, exec=ex, ms=1e3*ex/steps, pIter=it/steps, end=t)


def analyse(c):
    d = os.path.join(HERE, c)
    row = dict(case=c, **cost(d))
    if row["end"] > 29.99:
        m = F.metrics(F.load(HERE, c), REF, REF["t"][-1])
        row.update(dH=m["dHfit"], dc=m["dcMean"], **{k: v for k, v in F.fields(d).items()
                                                    if k in ("emax", "r0", "ulow", "Uair")})
    return row


if __name__ == "__main__":
    cases = sorted(os.path.basename(p) for p in glob.glob(HERE + "/*_B[0-3]"))
    with Pool(12) as p:
        rows = {r["case"]: r for r in p.map(analyse, cases)}
    keys = [("steps", "{:7d}"), ("exec", "{:8.0f}"), ("ms", "{:7.1f}"), ("pIter", "{:6.1f}"),
            ("dH", "{:+7.2f}"), ("dc", "{:+7.3f}"), ("emax", "{:7.2f}"), ("r0", "{:6.3f}"),
            ("ulow", "{:+7.4f}"), ("Uair", "{:6.2f}")]
    L = ["Pressure-solver study, max Co 0.1, 30 s. B0 papers' settings (PCG/DIC, nNonOrth 1,",
         "p_rghFinal 1e-12); B1 nNonOrth 0; B2 B1 + GAMG; B3 B2 + p_rghFinal 1e-9.",
         "exec: ExecutionTime [s]; ms: per step; pIter: pressure iterations per step;",
         "dH, dc, emax [%], r0, ulow, Uair as in papers/part2/finalSolver.py; 'x B0': exec/exec(B0).",
         "", f"{'case':24s}" + "".join(f"{k:>8s}" for k, _ in keys) + f"{'x B0':>7s}"]
    for c in cases:
        r, r0 = rows[c], rows[c[:-2] + "B0"]
        cells = []
        for k, fmt in keys:
            v = r.get(k)
            if v is not None and k == "emax":
                v = 100*v
            cells.append(f"{'-':>8s}" if v is None else fmt.format(v).rjust(8))
        L.append(f"{c:24s}" + "".join(cells) + f"{r['exec']/r0['exec']:7.2f}")
    open(os.path.join(HERE, "summary.txt"), "w").write("\n".join(L) + "\n")
    print("\n".join(L))
