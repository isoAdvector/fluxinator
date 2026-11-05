# geometricGravity

Code and test cases for the paper series *Consistent interfacial flow numerics*:
Part I, the gravity force term, and Part II, the momentum convection term. The
papers are the submodules `paper/` (Part I,
[isoAdvector/GravityFormulationPaper](https://github.com/isoAdvector/GravityFormulationPaper))
and `paper2/` (Part II,
[isoAdvector/interfaceConvectionPaper](https://github.com/isoAdvector/interfaceConvectionPaper)).

## The solver

interIsoGoam is interIsoFoam (isoAdvector, OpenFOAM v2512) with run-time options
for the gravity term (Part I) and the momentum convection term (Part II). It is
meant as a drop-in replacement for interIsoFoam:

| Settings | Behaves as |
|---|---|
| defaults, one outer corrector | interIsoFoam |
| `gravityModel geometric` (constant/gravityProperties) | the solver of Part I |
| also PIMPLE `explicitConvection yes` and `div(rhoPhi,U) Gauss laxWendroff U phi interfaceUpwind alpha.water 1e-6 limitedLinearV 1` | the final solver of Part II |
| `nOuterCorrectors` > 1 | time-centred step (not interIsoFoam, see below) |

The geometric gravity model takes the height in the gravity term at the
centroids of the reconstructed interface segments instead of at the face
centres. With more than one outer corrector the step is centred in time by
default: the interface is advected with the mean of the old and the latest
flux, and the convection term and the gravity force are the means of their
values at the start and the end of the step (PIMPLE `convectionTheta` and
`gravityTheta`, default 0.5; use `momentumPredictor no`). Stock interIsoFoam
centres only the interface advection, which damps waves. The header of
`applications/solvers/interIsoGoam/interIsoGoam.C` and `timeWeights.H` list
all switches; the diagnostic switch `momentumBudget` is off by default. Test
forces are set per case with fvOptions, not in the solver (for example the
uniform body acceleration of `run/periodicDisc`, a `vectorCodedSource` in
`system/fvOptions`).

## Contents

| Directory | Contents |
|---|---|
| `src/geometricVoF/` | isoAdvector and interface reconstruction used by interIsoGoam (library `newGeometricVoF`) |
| `src/gravityModels/` | gravity models: `basic` (as in OpenFOAM), `geometric` (Part I) and `constantHeight` (the exact height of a flat interface) |
| `src/interfaceSchemes/` | the interfaceUpwind and laxWendroff interpolation schemes (Part II); `libs (interfaceSchemes);` in controlDict |
| `applications/solvers/interIsoGoam/` | the solver |
| `applications/test/convectionTest/` | prescribed-flow test of the momentum convection (Part II) |
| `applications/test/ddtCorrSwitch/` | stock interFoam and interIsoFoam with a PIMPLE `ddtCorr` switch (Part I, Appendix A) |
| `applications/utilities/` | `setHydrostaticPressure`; `surfacePosition` (submodule) |
| `run/` | the test cases, one directory per case (`run/README.md`) |
| `papers/part1/`, `papers/part2/` | one directory per figure and table of each paper (below) |
| `paper/` | the Part I paper (submodule) |
| `paper2/` | the Part II paper (submodule) |

## Test cases, figures and tables

Each case in `run/` has a README, its base case(s), its studies (one directory
per study, named by what it varies, with the parameter file `<study>.csv` for
parmScanner studies) and `data/` with the summaries the papers cite.

Each figure and table of a paper has its own directory, named after its LaTeX
label: `papers/part1/fig_X/` for `fig:X` and `papers/part1/tab_X/` for `tab:X`
(likewise for Part II). It holds the script `make.py` (or `make.sh`; pvbatch for
ParaView figures) that makes the figure or table from the runs, and the file the
paper shows. The caption of each figure and table names its directory. Scripts
that several figures share are in `papers/part1/` and `papers/part2/`; a shared
script is run by the `make.py` of each of its outputs.

| Paper section | Case |
|---|---|
| I 2 The gravity term and its discrete balance | `run/progressiveWave` (nLayers study) |
| I 3 A flat interface at rest in a tilted box | `run/tiltedBox` |
| I 4 A flat interface at rest on non-orthogonal meshes | `run/skewedBox`, `run/refinedBox` |
| I 5 A standing wave across a refinement boundary | `run/standingWave` |
| I 6 A progressive wave | `run/progressiveWave` |
| I Appendix A The ddtCorr term | `<study>_ddtCorrOff` in every case |
| II 1, 6, 7 The progressive wave | `run/progressiveWave` (`convection/`, `finalSolver_*`) |
| II 3 The time lag | `run/taylorGreenVortex`, `run/progressiveWave` |
| II 4, 5 Prescribed flows, uniform acceleration | `run/prescribedFlows`, `run/periodicDisc` |

Shared Python modules (stream-function wave analysis, field reader, scan
summary) and helper scripts are in `run/common/`.

## Building

```shell
source <OpenFOAM-v2512>/etc/bashrc
./Allwmake
```

`applications/test/ddtCorrSwitch` is built separately (`wmake` in each solver
directory). The wave cases need `streamFunFoam` and `calcStreamFunctionWave`
from [InterFlowers/streamFunctionWave](https://github.com/InterFlowers/streamFunctionWave),
and the scans use `parmScanner` from
[isoAdvector/cfdTools](https://github.com/isoAdvector/cfdTools). The figure
scripts need python3 with numpy, scipy and matplotlib, and ParaView (pvbatch)
for the ParaView figures.

## Names before 2026-11-04

The repository was reorganised on 2026-11-04. Runs made before that date keep
the old names in their logs and in some generated files:

- The Part II options lived in a separate solver, `interIsoGoamExt`
  (`log.interIsoGoamExt`); it is merged into interIsoGoam, which gives bitwise
  identical results with the same settings. Exception: cases with
  `nOuterCorrectors` > 1 that relied on the default time weights ran with
  `convectionTheta 1` and `gravityTheta 1` unless `explicitConvection yes` was
  set. The PIMPLE entries `bodyAcceleration`, `bodyForceDensity` and
  `energyBudget` were removed (the disc cases now use `system/fvOptions`,
  bitwise identical).
- Directories and studies were moved and renamed by purpose; `run/README.md`
  lists the old names.

The full map of every move, the superseded code and the runs no paper uses are
outside this repository, in `../geometricGravity_archive` (`MOVES.txt`).
