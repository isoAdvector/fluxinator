# geometricGravity

Code and test cases for *Consistent interfacial flow numerics - Part I: The
gravity force term*. The paper itself is the submodule `paper/`
([isoAdvector/GravityFormulationPaper](https://github.com/isoAdvector/GravityFormulationPaper)).

The solver is interIsoGoam: interIsoFoam (isoAdvector, OpenFOAM v2512) with a
run-time selectable gravity model. The geometric gravity model takes the height
in the gravity term at the centroids of the reconstructed interface segments
instead of at the face centres.

## Contents

| Directory | Contents |
|---|---|
| `gravityModels/` | gravity models: `basic` (as in OpenFOAM), `geometric` (the paper's), `constantHeight`, `densityWeighted`, and the older variants `geometric.v3` and `gravityRecon` |
| `newGeometricVoF/` | isoAdvector and interface reconstruction used by interIsoGoam |
| `applications/interIsoGoam/` | the solver |
| `applications/ddtCorrSwitch/` | stock interFoam and interIsoFoam with a PIMPLE `ddtCorr` switch (Appendix A) |
| `applications/utilities/` | `setHydrostaticPressure`; `surfacePosition` (submodule) |
| `run/` | the test cases of the paper (below) |
| `archive/` | superseded code, kept for the record, not compiled |
| `paper/` | the paper (submodule) |

## Test cases and the paper

Every figure and table in the paper names the directory or file in this
repository its numbers come from. The case directories have a README and an
`Allrun` or scan script that rebuilds their results; `paperData/` holds the
summaries the paper cites. Generated cases and figures are not tracked.

| Paper section | Directory |
|---|---|
| 2 The gravity term and its discrete balance | `run/paper1/` (figure scripts), `run/periodicWaveFlume/` |
| 3 A flat interface at rest in a tilted box | `run/tiltedBox/` |
| 4 A flat interface at rest on non-orthogonal meshes | `run/skewBox/`, `run/refinedBox/` |
| 5 A standing wave across a refinement boundary | `run/standingWave/` |
| 6 A progressive wave | `run/periodicWaveFlume/`, `run/paper1/` |
| Appendix A The ddtCorr term | `run/ddtCorrOffPart1/` |

`run/qwistWave/` holds the stream-function wave analysis (`kinematics.py`,
`hmatch.py`) used by the wave figures, and the wave case of Qwist and
Christensen (2023). `run/summariseScan.py` and `run/snGradRhoReach` are shared
by the cases at rest.

## Building

```shell
source <OpenFOAM-v2512>/etc/bashrc
./Allwmake
```

`applications/ddtCorrSwitch` is built separately (`wmake` in each solver
directory). The wave cases need `streamFunFoam` and `calcStreamFunctionWave`
from [InterFlowers/streamFunctionWave](https://github.com/InterFlowers/streamFunctionWave),
and the scans use `parmScanner` from
[isoAdvector/cfdTools](https://github.com/isoAdvector/cfdTools).
