# Hydrostatic interface in a tilted box

A 1 m x 1 m two-dimensional box with 25 x 25 cells is rotated by 5 degrees about
its centre while gravity stays vertical. Water (1000 kg/m3) lies below the
horizontal plane z = 0.5 m and air above; both fluids are inviscid and surface
tension is zero. The exact solution is rest, U = 0, so every velocity is a
spurious current. Fixed time step 0.01 s, OpenFOAM v2512.

## Studies

| Study (directory, parameter file `<study>.csv`) | Content |
|---|---|
| `gravityModels_1step`, `_100steps`, `_1000steps` | interFoam (basic), interIsoFoam (interIsoGoam with basic), interIsoGoam with constantHeight and geometric, run for 1, 100 and 1000 time steps |
| `meshTimeStep_1step` | first time step only, dt = 1e-4 ... 1 s on 25^2, 50^2, 100^2 and 200^2 cells, interFoam and interIsoFoam |
| `airDensity_1000steps` | air density 1e-4 ... 1e3 kg/m3, 1000 time steps, interFoam, interIsoFoam and interIsoGoam (geometric) |
| `nLayers_1000steps` | interIsoGoam, geometric with nLayers 0, 1 and 2, and constantHeight, 1000 time steps |
| `<study>_ddtCorrOff` | the same runs with PIMPLE `ddtCorr no` (Part I, Appendix A); stock interFoam and interIsoFoam run as interFoamDdtCorr and interIsoFoamDdtCorr |

The base case `constant/gravityProperties` contains `nLayers 1;` so that the
studies can set it.

## Findings

- With the standard gravity term both solvers reach max|U| = 0.33 m/s after the
  first step, agreeing in all fifteen digits written and proportional to dt
  (32.97 dt m/s).
  After 10 s interFoam reaches 1.6 m/s and interIsoFoam 4.4 m/s.
- The constantHeight ground truth and the geometric model stay at 1e-10 to
  2e-9 m/s over 1000 steps.
- With lighter air the geometric model's residual grows by about one order of
  magnitude per decade of density, to 5.6e-5 m/s at 1e-4 kg/m3, while staying
  8 to 11 orders of magnitude below interIsoFoam. interFoam diverges at 1e-4 kg/m3.
- The spurious currents of the standard gravity term are concentrated in a few
  cells, and which cells those are follows in closed form from the tilt and the
  water level (closed-form theory in the archive, see below).
- The first ring of cells around the interface must get interface heights:
  with nLayers 0 (heights in interface cells only) the geometric model reaches
  0.53 m/s after the first step and 5.2 m/s after 10 s, worse than the standard
  gravity term. nLayers 1 and 2 stay below 2.5e-9 m/s, like constantHeight.

## Reproducing the results

```sh
./Allrun             # studies and data/
./Allclean [-scans]  # remove generated files (and with -scans the studies)
```

`writeData.py` writes a summary of each study to `data/` (tracked in git).
The paper figures and tables of this case are made in
`papers/part1/fig_tiltedBoxFields`, `fig_tiltedBoxRhoAir`, `tab_tiltedBox` and
`tab_firstStep`. `tiltedBoxCase.py` reads a case (fields, mesh vertices, tilt,
water level).

Analyses that no paper uses (the closed-form theory of which cells the standard
gravity term hits, the tilt and water-level scans, the interface band
statistics) are in `../../../geometricGravity_archive/reorg_20261104/run/tiltedBox`.
