# Test cases

| Case | Content | Papers |
|---|---|---|
| `tiltedBox` | flat interface at rest in a tilted box | Part I Sec. 3 |
| `skewedBox` | flat interface at rest on a skewed (non-orthogonal) mesh | Part I Sec. 4 |
| `refinedBox` | flat interface at rest in a box with a refined region | Part I Sec. 4 |
| `standingWave` | standing wave across a refinement boundary | Part I Sec. 5 |
| `progressiveWave` | stream function wave in a periodic domain of one wavelength | Part I Secs. 2, 6; Part II Secs. 1, 3, 6, 7 |
| `taylorGreenVortex` | translating Taylor-Green vortex, single phase (the time lag) | Part II Sec. 3 |
| `prescribedFlows` | water disc in air, prescribed velocity, momentum equation only | Part II Secs. 4, 5 |
| `periodicDisc` | water disc in a periodic box, full solver, uniform acceleration or no force | Part II Secs. 4, 5 |
| `common` | shared Python modules and helper scripts | |

Every case has a README. A study directory holds the runs of one parameter
scan or set of variants; `<study>_ddtCorrOff` repeats a Part I study with
PIMPLE `ddtCorr no` (Part I, Appendix A; `common/runStoredDdtCorrOff.sh`).
The figures and tables of the papers are made in `../papers`.

Runs are not tracked by git, only their set-ups and the summaries in `data/`;
the run data behind the papers is published separately.

## Old names (before 2026-11-04)

| Old | New |
|---|---|
| `run/periodicWaveFlume/gravityModelScan`, `gravityModelScan_fixedDt`, `nLayersScan` | `run/progressiveWave/solvers_maxCo`, `solvers_fixedDt`, `nLayers_maxCo` |
| `run/periodicWaveFlume/finalSolver/scan_maxCo`, `scan_fixedDt`, `IUonly` | `run/progressiveWave/finalSolver_maxCo`, `finalSolver_fixedDt`, `finalSolverNoLW` |
| `interfaceRelativePressure/run/waveLab/w*/...` | `run/progressiveWave/convection/<variant>_<mesh>_dt<dt>` |
| `run/ddtCorrOffPart1/<case>/<scan>` | `run/<case>/<study>_ddtCorrOff` |
| `run/tiltedBox/gravityModelScan_{1,100,1000}deltaT`, `meshDeltaTScan`, `airDensityScan_1000deltaT`, `nLayersScan` | `gravityModels_{1step,100steps,1000steps}`, `meshTimeStep_1step`, `airDensity_1000steps`, `nLayers_1000steps` |
| `run/skewBox/skewBoxScan_z`, `skewBoxScan`, `nonOrthCorrScan`, `hydrostaticInitScan` | `run/skewedBox/schemes_10steps`, `schemesRotatedMesh_10steps`, `nonOrthCorrectors_10steps`, `hydrostaticStart_10steps` |
| `run/refinedBox/scan` | `run/refinedBox/waterLevel_10steps` |
| `run/standingWave/scan` | `run/standingWave/solvers` |
| `interfaceRelativePressure/run/convectionTest/{uniform,accel,shearX,rot}*` | `run/prescribedFlows/{uniform,acceleration,shear,rotation}*` (`Evolve` is now `Carried`) |
| `interfaceRelativePressure/run/convectionTest/sv_*`, `u0_*`, `r4_*` | `run/periodicDisc/<convection>_rho{New,Old}`, `noForce_*`, `noForceRatio1e4_*` |
| `interfaceRelativePressure/run/tgvLag` | `run/taylorGreenVortex` |
| `paperData/` | `data/` |
| `run/paper1`, `interfaceRelativePressure/run/paper2` | `papers/part1`, `papers/part2` |

The complete list of moves is in
`../../geometricGravity_archive/reorg_20261104/moves.tsv`.
