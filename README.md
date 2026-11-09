# fluxinator\*

Developed and maintained by Johan Roenby, [STROMNING](https://stromning.com).

Code and test cases for the paper series *Consistent interfacial flow numerics*:

- J. Roenby, K. Missios, K. Møller, S. Aliyar and H. Scheufler (2026). Consistent
  interfacial flow numerics – Part I: The gravity force term.
  [arXiv:XXXX.XXXXX](https://arxiv.org/abs/XXXX.XXXXX)
- J. Roenby (2026). Consistent interfacial flow numerics – Part II: The momentum
  convection term. [arXiv:YYYY.YYYYY](https://arxiv.org/abs/YYYY.YYYYY)

Please cite the papers when you use the code.

STROMNING offers consulting, custom development, training and support for interfacial
and multiphase flow simulation with OpenFOAM. Contact: johan@stromning.com.

## The solver

**`fluxinator`** is a two-phase Volume-of-Fluid solver derived from OpenFOAM's
interIsoFoam. It keeps the geometric interface advection of interIsoFoam (isoAdvector)
and adds a consistent gravity force (Part I) and momentum convection (Part II). It is
built from the source of interIsoGoam (below), with other defaults
(`applications/solvers/interIsoGoam/solverDefaults.H`): the geometric gravity model, explicit convection
(PIMPLE `explicitConvection yes`) and the `ddtCorr` term off (PIMPLE `ddtCorr no`);
an entry in the case overrides them. The convection scheme is set in `system/fvSchemes`,
and explicit convection requires `interfaceUpwind` (the solver warns without it):

```
div(rhoPhi,U) Gauss laxWendroff interfaceUpwind limitedLinearV 1;
```

Keep the maximum Courant number at about 0.1 or less for accurate wave propagation.

`interIsoGoam` (`applications/solvers/interIsoGoam`) is interIsoFoam with the additions
of Parts I and II, all selected at run time; with the defaults of interIsoFoam and one outer
corrector it is interIsoFoam. It gives the same results as `fluxinator` with the same
settings. It is the solver of Part I and is kept so that the cases, which rely on its
defaults, reproduce the results of the papers. It has these run-time options:

| Settings | Behaves as |
|---|---|
| defaults, one outer corrector | interIsoFoam |
| `gravityModel geometric` (constant/gravityProperties) | the solver of Part I |
| also PIMPLE `explicitConvection yes` and `div(rhoPhi,U) Gauss laxWendroff interfaceUpwind limitedLinearV 1` | fluxinator of Part II |
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
all switches. Test forces and diagnostics are set per case, not in the solver:
forces with fvOptions (for example the uniform body acceleration of
`run/periodicDisc`, a `vectorCodedSource` in `system/fvOptions`), diagnostics
with function objects (for example `convectiveMomentumFlux` in
`run/progressiveWave/momentumFlux`).

## Contents

| Directory | Contents |
|---|---|
| `src/geometricVoF/` | isoAdvector and interface reconstruction: a modified copy of OpenFOAM's geometricVoF, in namespace `Foam::isoAdvection` so that both can be loaded together (library `newGeometricVoF`; sampled surface `reconstructedInterface`) |
| `src/gravityModels/` | gravity models: `basic` (as in OpenFOAM), `geometric` (Part I) and `constantHeight` (the exact height of a flat interface) |
| `src/interfaceSchemes/` | the interfaceUpwind and laxWendroff interpolation schemes (Part II); fluxinator links the library, other solvers need `libs (interfaceSchemes);` in controlDict |
| `src/functionObjects/` | the function object `convectiveMomentumFlux`, the face momentum flux of the convection term (Part II); `libs (interfaceFunctionObjects);` |
| `applications/solvers/fluxinator/` | fluxinator: the source of interIsoGoam with the defaults of `solverDefaults.H` |
| `applications/solvers/interIsoGoam/` | interIsoGoam, and the source of both solvers |
| `applications/test/convectionTest/` | prescribed-flow test of the momentum convection (Part II) |
| `applications/test/ddtCorrSwitch/` | stock interFoam and interIsoFoam with a PIMPLE `ddtCorr` switch (Part I, Appendix A) |
| `applications/utilities/` | `setHydrostaticPressure` (hydrostatic `p_rgh` of a flat interface at rest); `surfacePosition` (submodule; function object for the interface height, a geometric alternative to `interfaceHeight`) |
| `run/` | the test cases, one directory per case (`run/README.md`) |
| `papers/part1/`, `papers/part2/` | one directory per figure and table of each paper (below) |

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
source <OpenFOAM-v2606>/etc/bashrc
./Allwmake
```

After cloning, get `surfacePosition` with
`git submodule update --init applications/utilities/surfacePosition`. The wave cases
need `streamFunFoam` and `calcStreamFunctionWave` from [InterFlowers/streamFunctionWave](https://github.com/InterFlowers/streamFunctionWave),
and the scans use `parmScanner` from
[isoAdvector/cfdTools](https://github.com/isoAdvector/cfdTools). The figure
scripts need python3 with numpy, scipy and matplotlib, and ParaView (pvbatch)
for the ParaView figures.

## Running

`run/quickTest/Allrun` builds nothing and takes about a minute: it runs
fluxinator on the accelerated disc and on one period of the wave, and prints the results
next to the expected values. Each case in `run/` has an `Allrun` for its studies
(`run/README.md`); `MAXJOBS` limits the number of simultaneous runs. All scripts use
the OpenFOAM environment of the shell they are started from.

## Names in old logs

Logs of some Part II runs name the solver `interIsoGoamExt`, an earlier version of
interIsoGoam that gives bitwise identical results with the same settings. The
`budget_*` fields in `run/progressiveWave/convection` were written by a diagnostic
that is now the function object `convectiveMomentumFlux`.

## Use of generative AI

The author used Claude Code with Anthropic's Claude Opus models extensively in the
development of this repository, including the solver, libraries and utilities, the case
setups, and the scripts for running, post-processing and making the figures and tables.
All code, cases and scripts produced with the tool have been thoroughly reviewed by the
author, who takes full responsibility for them.

## Licence

The code is released under the GNU General Public License v3 (`LICENSE`), as OpenFOAM.

---

\* Not affiliated with or endorsed by Dr. Heinz Doofenshmirtz or Doofenshmirtz Evil
Incorporated.
