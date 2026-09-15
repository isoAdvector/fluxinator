# Hydrostatic interface in a tilted box

A 1 m x 1 m two-dimensional box with 25 x 25 cells is rotated by 5 degrees about
its centre while gravity stays vertical. Water (1000 kg/m3) lies below the
horizontal plane z = 0.5 m and air above; both fluids are inviscid and surface
tension is zero. The exact solution is rest, U = 0, so every velocity is a
spurious current. Fixed time step 0.01 s, OpenFOAM v2512.

## Scans

| Parameter file | Content |
|---|---|
| `gravityModelScan_1deltaT.csv`, `_100deltaT.csv`, `_1000deltaT.csv` | interFoam (basic), interIsoFoam (interIsoGoam with basic), interIsoGoam with constantHeight, geometric and gravityRecon, run for 1, 100 and 1000 time steps |
| `deltaTScan.csv` | first time step only, dt = 1e-4 ... 1 s, interFoam and interIsoFoam |
| `airDensityScan_1000deltaT.csv` | air density 1e-4 ... 1e3 kg/m3, 1000 time steps, interFoam, interIsoFoam and interIsoGoam (geometric) |

## Findings

- With the standard gravity term both solvers reach max|U| = 0.33 m/s after the
  first step, identical to seven digits and proportional to dt (32.97 dt m/s).
  After 10 s interFoam reaches 1.6 m/s and interIsoFoam 4.4 m/s.
- The constantHeight ground truth and the geometric model stay at 1e-10 to
  2e-9 m/s over 1000 steps.
- With lighter air the geometric model's residual grows by about one order of
  magnitude per decade of density, to 5.6e-5 m/s at 1e-4 kg/m3, while staying
  8 to 11 orders of magnitude below interIsoFoam. interFoam diverges at 1e-4 kg/m3.

## Reproducing the figures

```sh
./Allrun             # scans and figures
./Allclean [-scans]  # remove generated files (and with -scans the scans)
```

Figures: `figure_paper.pdf` (velocity fields, from `figure_paper.tex`),
`maxU_tiltedBox.pdf` (`plotMaxU.py`), `firstStepMaxU_tiltedBox.pdf`
(`plotFirstStepMaxU.py`) and `maxUvsRhoAir_tiltedBox.pdf` (`plotMaxUvsRhoAir.py`).

The interFoam logs end with a fatal error from the interface sampling function
object in `controlDict`, which needs a reconstruction scheme that interFoam does
not have. It occurs after the last time step and does not affect the results.
