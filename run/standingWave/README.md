# Standing wave on a locally refined mesh

A linear standing wave (k = pi, depth 0.5 m) in a 2-D tank, on a uniform mesh and
on a mesh refined around the interface, with interFoam, interIsoFoam and
interIsoGoam (geometric gravity model, nLayers 2). The wave starts at rest at
maximum elevation, so the exact velocity after one step of 0.01 s is only
0.0154 m/s and anything above that is error. Part I, Section 5.

| Directory | Content |
|---|---|
| `baseCase` | the case set-up (`Allrun` builds the refined mesh and runs it) |
| `baseUniformMesh`, `baseRefinedMesh` | prepared base cases copied by `runCase.sh`; `makeBaseMeshes.sh` builds their meshes |
| `solvers` | `paperScan.sh`: every solver on both meshes at twelve interface offsets over one coarse cell (`scan.csv`), and two-period histories on the refined mesh (`history_<solver>.csv`, cases `hist_<solver>`) |
| `snapshots` | `snapshots.sh`: the three solvers on the refined mesh with dt = T/118, written at T/2 and T |
| `<study>_ddtCorrOff` | the same runs with PIMPLE `ddtCorr no` (Part I, Appendix A; `*_ddtCorrOff.sh`) |

`runCase.sh <outdir> [key=value ...]` copies a base case (`$SWBASE`, default
`baseRefinedMesh`), applies the settings and runs it.

The paper figure and table are made in `papers/part1/fig_standingWaveSnapshots`
and `papers/part1/tab_standingWave`. The experiments that located the remaining
spurious currents (`FINDINGS.md`) are not part of the release.
