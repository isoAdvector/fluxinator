# Hydrostatic interface in a skewed box

A 1 m x 1 m box with a parallelogram mesh of 25 x 25 cells: the cell rows slope
by 0.1 m over the 1 m width (5.7 degrees non-orthogonality). The mesh is rotated
about the box centre by `rotAngle` in `rotAngleDict`: 90 degrees in `baseCase`
(rows level, columns sheared), and the studies except `schemesRotatedMesh_10steps` set
`rotAngle 0` (rows sloping relative to the horizontal interface). Water below z = 0.5 m, air above,
inviscid, no surface tension, dt = 0.01 s, 10 time steps, plicRDF, no
non-orthogonal corrector loops unless stated. Exact solution U = 0.

`baseCaseHydrostatic` is identical except that the `setHydrostaticPressure`
utility (`../../applications/utilities/setHydrostaticPressure`) sets a hydrostatic initial
p_rgh = -rho (g & x_I) after `setAlphaField`. Both base
cases contain `nLayers 1;` in `constant/gravityProperties` so that the studies can
set it.

## Studies

| Study (directory, parameter file `<study>.csv`) | Content |
|---|---|
| `schemes_10steps` | paper figure: interFoam, interIsoFoam (interIsoGoam with basic) and interIsoGoam (geometric), with corrected and uncorrected laplacian and snGrad schemes |
| `schemesRotatedMesh_10steps` | the same with the mesh rotated by 90 degrees (rows level) |
| `nonOrthCorrectors_10steps` | corrected schemes, nNonOrthogonalCorrectors 0, 1, 2, 5, 10 (geometric), 0, 5 (basic) and 0, 2 (constantHeight) |
| `hydrostaticStart_10steps` | hydrostatic initial p_rgh: geometric and constantHeight with 0 and 5 correctors, basic with 0 |
| `<study>_ddtCorrOff` | the same runs with PIMPLE `ddtCorr no` (Part I, Appendix A) |

## Findings (max|U| after 10 steps)

| | corrected | uncorrected |
|---|---|---|
| interFoam, basic | 5.6 m/s | 0.76 m/s |
| interIsoFoam, basic | 7.8 m/s | 1.3 m/s |
| interIsoGoam, geometric | 3.5e-2 m/s | 4.1e-11 m/s |

- With level rows (`schemesRotatedMesh_10steps`) uncorrected schemes alone keep all solvers at
  rest, because the faces carrying the density jump share one height.
- The residual of the geometric model with corrected schemes is a start-up
  error: p_rgh starts at 0 and the explicit non-orthogonal correction lags. Two
  non-orthogonal correctors reduce it to 2e-8 m/s, and a hydrostatic initial
  p_rgh removes it (7.7e-10 m/s) without correctors. constantHeight, with the
  exact interface height, has the same start-up error (4.3e-2 m/s after the first
  step, 2.1e-8 m/s with 2 correctors).
- The diagnostic `../common/snGradRhoReach` shows why nLayers 1 suffices here: the
  corrected snGrad(rho) is non-zero on 128 faces (76 uncorrected), but all their
  cells lie in rings 0 and 1, because opposite cell faces are parallel and the
  correction vector is orthogonal to the face normal.
- nLayers 2 gives the same results as nLayers 1 on this mesh (scans in the
  archive, `geometricGravity_archive/reorg_20261104/run/skewedBox`).

## Running

```sh
./Allrun             # studies, case images and data/ (tracked in git)
./Allclean [-scans]  # remove generated files (and with -scans the studies)
```

The paper figure and table of this case are made in `papers/part1/fig_skewedBox`
and `papers/part1/tab_nonOrthSkew`.
