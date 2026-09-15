# Hydrostatic interface in a tilted box with a refined region

The tilted box of `../tiltedBox` (25 x 25 cells, rotated by 5 degrees) with the
square 0.25 < x, z < 0.75 refined once (2 x 2 in x and z, `topoSetDict1`,
`refineMeshDict1`) before the rotation. The refined region straddles the water
line. The mesh has 1132 cells, including 52 polyhedra at the refinement
boundary, with non-orthogonality up to 18.4 degrees and skewness 0.33.

`baseCaseHydrostatic` is identical except that `setHydrostaticPrgh.py` sets a
hydrostatic initial p_rgh = -rho (g & x_I) after `setAlphaField`. Both base
cases contain `nLayers 1;` in `constant/gravityProperties` so that the scans can
set it. 10 time steps of 0.01 s, corrected schemes unless stated.

## Scans

| Parameter file | Content |
|---|---|
| `zeroStartScan.csv` | p_rgh = 0 initially: geometric with 0, 1, 2, 5 non-orthogonal correctors, constantHeight with 0, 5, basic with 0, geometric with uncorrected schemes, constantHeight with 2 |
| `hydrostaticInitScan.csv` | hydrostatic start: geometric and constantHeight with 0, 5 correctors, basic with 0 |
| `nLayersZeroStartScan.csv` | geometric, nLayers 1 (0 correctors) and nLayers 2 (2 and 0 correctors) |
| `nLayersHydrostaticScan.csv` | geometric, nLayers 1, 2 and 3, hydrostatic start |

## Findings (max|U| after 10 steps)

| | zero start, 0 correctors | zero start, 2 or 5 correctors | hydrostatic start |
|---|---|---|---|
| basic | 3.3 m/s | | 3.3 m/s |
| constantHeight | 0.14 m/s | 1.4e-7 m/s (2), 1.1e-7 m/s (5) | 2.4e-9 m/s |
| geometric, nLayers 1 | 0.95 m/s | 0.95 m/s (2, 5) | 0.95 m/s |
| geometric, nLayers 2 | 0.15 m/s | 5.8e-8 m/s (2) | 6.4e-9 m/s |

- With corrected schemes, fvc::snGrad(rho) reaches one ring of cells further out
  than with uncorrected schemes. With nLayers 1 the geometric model does not
  give these faces the interface height, and the hydrostatic state is not a
  discrete equilibrium (first p_rgh residual 1.5e-5 from a hydrostatic start).
- With nLayers 2 the geometric model behaves like the exact constantHeight
  model (first residual 3e-13); nLayers 3 gives no further improvement. Like
  constantHeight it still needs a hydrostatic start or non-orthogonal correctors:
  from a zero start without correctors it reaches 0.15 m/s.
- With uncorrected schemes nLayers 1 suffices (1.3e-9 m/s).
- The diagnostic `../snGradRhoReach` (a coded function object for
  `postProcess -time 0 -dict`) marks the cells adjacent to faces with a non-zero
  snGrad(rho) and prints their ring levels. Uncorrected: 122 faces, all cells in
  rings 0 and 1. Corrected: 130 faces, two cells in ring 2 where the interface
  crosses the refinement boundary. The maximum velocity with nLayers 1 (0.95 m/s)
  is in one of these two cells.

## Running

```sh
./Allrun             # scans, summary_*.txt and paper figures
./Allclean [-scans]  # remove generated files (and with -scans the scans)
```

Allrun also renders `|U|` images of the hydrostatic start cases (basic,
constantHeight, geometric with nLayers 1 and 2) with `generateFigure.py` and
collects them in `figure_tiltedBoxRefined.pdf`
(`figure_tiltedBoxRefined.tex`), and renders `paperFigures/reach_uncorrected.png` and `reach_corrected.png` with `renderSnGradReach.py` (pvbatch), and
plots `../nonOrthMaxU.pdf` with `../plotNonOrthMaxU.py`, which also needs the scans
in `../skewBox` (run `../skewBox/Allrun` first).
