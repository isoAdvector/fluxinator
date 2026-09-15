# Hydrostatic interface in a skewed box

A 1 m x 1 m box with a parallelogram mesh of 25 x 25 cells: the cell rows slope
by 0.1 m over the 1 m width (5.7 degrees non-orthogonality). The mesh is rotated
about the box centre by `rotAngle` in `rotAngleDict`: 90 degrees in `baseCase`
(rows level, columns sheared), and the scans except `skewBoxScan` set
`rotAngle 0` (rows sloping relative to the horizontal interface). Water below z = 0.5 m, air above,
inviscid, no surface tension, dt = 0.01 s, 10 time steps, plicRDF, no
non-orthogonal corrector loops unless stated. Exact solution U = 0.

`baseCaseHydrostatic` is identical except that `setHydrostaticPrgh.py` sets a
hydrostatic initial p_rgh = -rho (g & x_I) after `setAlphaField`.

## Scans

| Parameter file | Content |
|---|---|
| `skewBoxScan_z.csv` | paper figure: interFoam, interIsoFoam (interIsoGoam with basic) and interIsoGoam (geometric), with corrected and uncorrected laplacian and snGrad schemes |
| `skewBoxScan.csv` | the same with the mesh rotated by 90 degrees (rows level) |
| `nonOrthCorrScan.csv` | corrected schemes, nNonOrthogonalCorrectors 0, 1, 2, 5, 10 (geometric) and 0, 5 (basic) |
| `hydrostaticInitScan.csv` | hydrostatic initial p_rgh: geometric and constantHeight with 0 and 5 correctors, basic with 0 |
| `nLayersDefaultScan.csv` | geometric without an nLayers entry (default), 2 correctors |

## Findings (max|U| after 10 steps)

| | corrected | uncorrected |
|---|---|---|
| interFoam, basic | 5.6 m/s | 0.76 m/s |
| interIsoFoam, basic | 7.8 m/s | 1.3 m/s |
| interIsoGoam, geometric | 3.5e-2 m/s | 4.1e-11 m/s |

- With level rows (`skewBoxScan`) uncorrected schemes alone keep all solvers at
  rest, because the faces carrying the density jump share one height.
- The residual of the geometric model with corrected schemes is a start-up
  error: p_rgh starts at 0 and the explicit non-orthogonal correction lags. Two
  non-orthogonal correctors reduce it to 2e-8 m/s, and a hydrostatic initial
  p_rgh removes it (7.7e-10 m/s) without correctors.

## Running

```sh
./Allrun             # scans, figure_skewBox.pdf and summary_*.txt
./Allclean [-scans]  # remove generated files (and with -scans the scans)
```
