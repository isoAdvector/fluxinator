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
| `meshDeltaTScan.csv` | first time step only, dt = 1e-4 ... 1 s on 25^2, 50^2, 100^2 and 200^2 cells, interFoam and interIsoFoam |
| `airDensityScan_1000deltaT.csv` | air density 1e-4 ... 1e3 kg/m3, 1000 time steps, interFoam, interIsoFoam and interIsoGoam (geometric) |
| `nLayersScan.csv` | interIsoGoam, geometric with nLayers 0, 1 and 2, and constantHeight, 1000 time steps |
| `tiltScan_1deltaT.csv` | interFoam, first time step only, box tilt 2 to 20 degrees in steps of 0.1 |
| `waterLevelScan_1deltaT.csv` | interFoam, first time step only, water level +/- 0.045 m at tilt 5, 10 and 15 degrees |

The base case `constant/gravityProperties` contains `nLayers 1;` so that the
scans can set it.

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
  water level: see `vertexOffsetTheory.py` and the summary below.
- The first ring of cells around the interface must get interface heights:
  with nLayers 0 (heights in interface cells only) the geometric model reaches
  0.53 m/s after the first step and 5.2 m/s after 10 s, worse than the standard
  gravity term. nLayers 1 and 2 stay below 2.5e-9 m/s, like constantHeight.

## Reproducing the figures

```sh
./Allrun             # scans and figures
./Allclean [-scans]  # remove generated files (and with -scans the scans)
```

The numbers behind the paper are written to `paperData/` (tracked in git) by
`paperData.py`: scan summaries, the max|U| table, the first step against dt and
the interface band and numerical vapour statistics (`analyseInterfaceBand.py`).

Figures: `figure_paper.pdf` (velocity fields, from `figure_paper.tex`),
`maxU_tiltedBox.pdf` (`plotMaxU.py`), `firstStepMaxU_tiltedBox.pdf`
(`plotFirstStepMaxU.py`), `maxUvsRhoAir_tiltedBox.pdf` (`plotMaxUvsRhoAir.py`)
`vertexOffset_tiltedBox.pdf` (`plotVertexOffset.py`) and
`offsetGeometry_tiltedBox.pdf` (`plotOffsetGeometry.py`, a definition sketch of
the symbols, pure geometry).

## Which cells are hit

`vertexOffsetTheory.py` holds a closed-form theory for the first time step with
the standard gravity term, and checks it against the interFoam runs of the two
scans above. It needs no simulation of its own.

Starting from rest the pressure equation removes exactly the discrete gradient
part of the face gravity force, so the velocity is driven only by what is left:
the imbalance of that force around the four-cell loops, one scalar per mesh
vertex. For a flat interface in a tilted square mesh that loop imbalance is a
universal function of the single offset `xi = z_w - z_vertex`: it vanishes beyond
one cell from the interface, peaks when the interface passes exactly through the
vertex, and the vertex offsets form a lattice, so the hot columns are those where
`I tan(theta)` comes closest to an integer and recur every `cot(theta)` columns.
Because the velocity scales as `dt/rho`, the worst vertex is not the one on the
interface but the one just on the light-fluid side, whose four cells are all
light: that window has width `tan(theta)` in fractional offset, and at its edge the
imbalance is still half its peak.

```sh
./vertexOffsetTheory.py            # identities, then the 238 interFoam runs
```

Over those runs the light-fluid window selects between none and nine vertices
per case, in line with `N tan(theta)`, and the cell carrying the largest
spurious velocity touches one of them in 94.5% of the runs. Ranking them singles
out one vertex, whose four cells contain the peak cell in 78% of the runs; the
predicted magnitude has a median ratio of 1.03 to the measured one. The output
is written to `paperData/vertexOffsetTheory.txt` by `paperData.py`.

`tiltedBoxCase.py` reads a case (fields, mesh vertices, tilt, water level) and
is used by `vertexOffsetTheory.py` and `plotVertexOffset.py`.

`firstStepTwin.py` is a numpy re-implementation of the first time step that
reproduces OpenFOAM to about 1e-10. No paper result depends on it; it is kept
because it lets the first step of the argument be checked directly, by replacing
the face force with a different one that has the same loop imbalance - something
stock OpenFOAM cannot do.
