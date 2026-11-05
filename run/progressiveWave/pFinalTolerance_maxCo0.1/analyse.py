"""p_rghFinal tolerance study: metrics as papers/part2/finalSolver.py, solver log data, and the
wave-gauge difference to the 1e-12 run. Writes summary.txt. Run with python3 analyse.py."""
import glob, os, re, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "../../../papers/part2"))
import finalSolver as F  # noqa: E402

REF = F.load(os.path.join(F.WAVE, "solvers_fixedDt"), "0036")
P = os.path.join(HERE, "../pressureSolver_maxCo0.1")
CASES = [("1e-9 (GAMG, nNonOrth 0)", P, "final_190x40_B3"), ("1e-10", HERE, "final_190x40_tol1e-10"),
         ("1e-11", HERE, "final_190x40_tol1e-11"), ("1e-12", P, "final_190x40_B0"),
         ("1e-13", HERE, "final_190x40_tol1e-13")]


def height(base, c):
    return np.loadtxt(os.path.join(base, c, "postProcessing/interfaceHeight1/0/height.dat"),
                      comments="#", usecols=(0, 2), unpack=True)


t0, h0 = height(P, "final_190x40_B0")
L = ["fluxinator, 190x40, max Co 0.1, 30 s; surfCellTol 1e-10. dH, dc, emax [%], r0 as finalSolver.py;",
     "maxDh: max |h - h(1e-12)| at the wave gauge [m]; contErr: mean of 'sum local' at the final corrector;",
     "dV: change of the water volume fraction from first to last step; warn: cutFaceAdvect warnings.", "",
     f"{'p_rghFinal':24s}{'steps':>7s}{'dH':>8s}{'dc':>8s}{'emax':>7s}{'r0':>7s}{'maxDh':>10s}{'contErr':>10s}{'dV':>10s}{'warn':>7s}"]
for name, base, c in CASES:
    m = F.metrics(F.load(base, c), REF, REF["t"][-1])
    f = F.fields(os.path.join(base, c))
    t, h = height(base, c)
    dh = np.abs(np.interp(t0, t, h) - h0).max()
    log = glob.glob(os.path.join(base, c, "log.interIsoGoam"))[0]
    steps, warn, vol, cont, last = 0, 0, [], [], None
    for line in open(log):
        if line.startswith("Time = "):
            steps += 1
            if last is not None:
                cont.append(last)
        elif "FOAM Warning" in line:
            warn += 1
        elif line.startswith("Phase-1 volume fraction"):
            vol.append(float(line.split()[4]))
        elif line.startswith("time step continuity errors"):
            last = float(re.search(r"sum local = (\S+),", line)[1])
    L.append(f"{name:24s}{steps:7d}{m['dHfit']:+8.2f}{m['dcMean']:+8.3f}{100*f['emax']:7.2f}{f['r0']:7.3f}"
             f"{dh:10.2e}{np.mean(cont):10.2e}{vol[-1]-vol[0]:+10.1e}{warn:7d}")
open(os.path.join(HERE, "summary.txt"), "w").write("\n".join(L) + "\n")
print("\n".join(L))
