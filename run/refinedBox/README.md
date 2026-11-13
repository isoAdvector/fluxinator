# Hydrostatic interface in a box with a refined region

The box of `../tiltedBox` (25 x 25 cells), not tilted, with the central 13 x 13 cells
(0.24 < x, z < 0.76) split 2 x 2 (`topoSetDict1`, `refineMeshDict1`). The mesh is
orthogonal except at the refinement boundary, where each coarse face towards the
refined region is split in two and is non-orthogonal by atan(1/3) = 18.4 degrees.
1132 cells, including 52 polyhedra, skewness 0.33.

The still water level z_w is scanned over one coarse cell, 0.48 ... 0.52 in steps of
0.0025 (17 levels; coarse faces at 0.48 and 0.52, fine faces also at 0.50), for the
configurations of `makeCases.sh`: interFoam and interIsoFoam (interIsoGoam, basic),
constantHeight, and the geometric model with nLayers 1, 2 (and 3), each with corrected
schemes from a zero p_rgh with 0 or 2 non-orthogonal correctors, from a hydrostatic
p_rgh (`setHydrostaticPressure`), and with uncorrected schemes. 10 time steps of 0.01 s.
Case names: `waterLevel_10steps/<configuration>_z<z_w>`.

In the interFoam cases the `surfaces` (interface) function object is removed: it needs
the isoAdvector reconstruction and aborts interFoam at the first write.

## Findings

- Standard gravity term: 1-6 m/s after one step with corrected schemes at every level.
  Uncorrected: round-off only when z_w lies on coarse faces; otherwise 0.6-4 m/s for
  levels below the centre of the cut coarse cells and 0.004-0.1 m/s above. The maximum
  is in the empty fine cell at the refinement boundary just above the interface: the
  unbalanced split faces border air.
- constantHeight: start-up error (0.01-0.5 m/s) from a zero start with corrected
  schemes, reduced to 1e-8 - 3e-7 m/s by 2 correctors and removed by the hydrostatic start
  (round-off, first residual 1e-14); round-off with uncorrected schemes. At z_w = 0.49, 5 and
  10 correctors give the same as 2, also for nLayers 2 (`data/nonOrthCorrectors_z0.49.txt`,
  cases `{cH,g2}_nc{1,5,10}_z0.4900`, run by `Allrun` after the scan).
- geometric, nLayers 1, corrected: 0.2-1 m/s at every level where the interface cuts
  cells, hydrostatic start or not (first residual 6e-5 - 2e-4); maximum in one of the
  four layer-2 cells of the corrected snGrad(rho) reach. nLayers 2 behaves like
  constantHeight at every level with the hydrostatic start; from a zero start with 2
  correctors it reaches 3.6e-5 m/s at z_w = 0.50 (constantHeight 2.9e-7); nLayers 3 gives the same as 2. Uncorrected: nLayers 1
  suffices.
- Interface on coarse cell faces (z_w = 0.48, 0.52): volume fractions are exactly 0 and
  1 and there are no interface cells. (At z_w = 0.50 the interface lies on fine faces in
  the refined region but cuts the coarse cells outside it in half.) The geometric model assigns the face height to both
  cells of every face between a full and an empty cell; without that rule it reduces to
  the basic model there (4-5 m/s with corrected schemes, any nLayers). From a zero start
  with 2 correctors the first-step velocity (2e-8 m/s) later perturbs some volume
  fractions by ~1e-9, creating nearly empty interface cells, and the velocity reaches
  1e-4 m/s after 10 steps; from the hydrostatic start it stays at round-off.

## Running

```sh
./Allrun             # study waterLevel_10steps (357 cases, about a minute), data/ and case images
./Allclean [-scans]  # remove generated files (and with -scans the study)
```

`Allrun` writes `data/waterLevel_10steps.txt` (`summarise.py`: first- and tenth-step
max|U| per configuration and level), `snGradRhoReach.txt`, `maxUCell.txt` (`maxUCell.py`)
and `checkMesh.txt`, and renders the case images of the velocity-field figure
(`generateFigure.py` in four cases at z_w = 0.49, pvbatch). The paper figures and the
table are made in `papers/part1/fig_reach`, `fig_refinedBox` and `tab_nonOrthRefined`
(at z_w = 0.49, with the range over all levels at which the interface cuts cells). `waterLevel_10steps_ddtCorrOff` holds the same runs with PIMPLE
`ddtCorr no` (Part I, Appendix A; `runCase_ddtCorrOff.sh`).
