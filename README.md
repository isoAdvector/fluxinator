# fluxinator

fluxinator\* is an incompressible interfacial flow solver for OpenFOAM based on the
isoAdvector geometric Volume-of-Fluid method. It is essentially the interIsoFoam solver,
but with numerics that take the sharp jump in density at the fluid interface into account.

Code and test cases for the paper series *Consistent interfacial flow numerics*:

- J. Roenby, K. Missios, K. Møller, S. Aliyar and H. Scheufler (2026). Consistent
  interfacial flow numerics – Part I: The gravity force term. arXiv link to follow.
- J. Roenby (2026). Consistent interfacial flow numerics – Part II: The momentum
  convection term. arXiv link to follow.

Please cite the papers when you use the code.

The code is actively developed and maintained by Johan Roenby, [STROMNING](https://stromning.com).
STROMNING offers consulting, custom development, training and support for interfacial
and multiphase flow simulation with OpenFOAM. Contact: johan@stromning.com.

## The solvers

The code has two solvers, built from the same source, that differ only in their
defaults:

| Solver | Defaults | Paper |
|---|---|---|
| **`fluxinator`** | the geometric gravity force and explicit momentum convection: the recommended settings | Part II |
| `interIsoGoam` | those of OpenFOAM's interIsoFoam; the additions are switched on in the case | Part I, and the steps of Part II |

Use `fluxinator`. `interIsoGoam` is maintained alongside it: with the geometric gravity
model it is the solver of Part I, and the cases of both papers use it to switch the
additions on one at a time.

**`fluxinator`** is a two-phase Volume-of-Fluid solver derived from OpenFOAM's
interIsoFoam. It keeps the geometric interface advection of interIsoFoam (isoAdvector)
and adds a consistent gravity force (Part I) and momentum convection (Part II). Its
source is in `applications/solvers/fluxinator`, and its defaults
(`solverDefaults.H` there) are the geometric gravity model, explicit convection
(PIMPLE `convectionNewTimeWeight 0`), the `ddtCorr` term off (PIMPLE `ddtCorr no`) and
three reconstructions of the initial interface (`nInitialReconstructions 3` in the solver
dictionary of alpha); an entry in the case overrides them. The convection scheme is set in `system/fvSchemes`,
and explicit convection requires `interfaceUpwind` (the solver warns without it):

```
div(rhoPhi,U) Gauss laxWendroff interfaceUpwind limitedLinearV 1;
```

Keep the maximum Courant number at about 0.1 or less for accurate wave propagation.

`interIsoGoam` (`applications/solvers/interIsoGoam`) is built from the source of fluxinator
with the defaults of interIsoFoam, so the additions of Parts I and II are all selected at run
time; with these defaults and one outer corrector it is interIsoFoam. It gives the same results as `fluxinator` with the same
settings. The cases of the papers rely on its defaults, so they reproduce the published
results unchanged. It has these run-time options:

| Settings | Behaves as |
|---|---|
| defaults, one outer corrector | interIsoFoam |
| `interfaceForces { gravity { model geometric; } }` (system/fvSchemes) | the solver of Part I |
| the defaults of fluxinator (above) and `div(rhoPhi,U) Gauss laxWendroff interfaceUpwind limitedLinearV 1` | fluxinator of Part II |
| `nOuterCorrectors` > 1 | time-centred step (not interIsoFoam, see below) |

The geometric gravity model takes the height in the gravity term at the
centroids of the reconstructed interface segments instead of at the face
centres. With more than one outer corrector the step is centred in time by
default: the interface is advected with the mean of the old and the latest
velocity and flux, and the convection term and the gravity force are the means of
their values at the start and the end of the step. The three weights are
`convectionNewTimeWeight` and `gravityNewTimeWeight` (PIMPLE) and `outerNewTimeWeight`
(solver dictionary of alpha), each the weight of the end of the step, default 0.5 with
outer correctors. They belong together: changing one of them alone can damp a wave
strongly. Stock interIsoFoam centres only the interface advection, which damps waves. The header of
`applications/solvers/fluxinator/fluxinator.C` and `readControls.H` list
all switches. Test forces and diagnostics are set per case, not in the solver:
forces with fvOptions (for example the uniform body acceleration of
`run/periodicDisc`, a `vectorCodedSource` in `system/fvOptions`), diagnostics
with function objects (for example `convectiveMomentumFlux` in
`run/progressiveWave/momentumFlux`).

## Contents

| Directory | Contents |
|---|---|
| `src/fluxinatorVoF/` | isoAdvector and interface reconstruction: a modified copy of OpenFOAM's geometricVoF, in namespace `Foam::fluxinatorVoF` so that both can be loaded together (library `fluxinatorVoF`; sampled surface `reconstructedInterface`) |
| `src/gravityModels/` | gravity models: `basic` (as in OpenFOAM), `geometric` (Part I) and `constantHeight` (the exact height of a flat interface) |
| `src/interfaceSchemes/` | the interfaceUpwind and laxWendroff interpolation schemes (Part II); fluxinator links the library, other solvers need `libs (interfaceSchemes);` in controlDict |
| `src/functionObjects/` | the function object `convectiveMomentumFlux`, the face momentum flux of the convection term (Part II); `libs (interfaceFunctionObjects);` |
| `applications/solvers/fluxinator/` | fluxinator, and the source of both solvers |
| `applications/solvers/interIsoGoam/` | interIsoGoam: the source of fluxinator with the defaults of interIsoFoam (`solverDefaults.H`) |
| `applications/test/convectionTest/` | prescribed-flow test of the momentum convection (Part II) |
| `applications/test/ddtCorrSwitch/` | stock interFoam and interIsoFoam with a PIMPLE `ddtCorr` switch (Part I, Appendix A) |
| `applications/utilities/` | `setHydrostaticPressure` (hydrostatic `p_rgh` of a flat interface at rest) |
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

The code is written for OpenFOAM.com v2606 and builds with v2306 and later. With
any version, edit the code in `src` and `applications` and build with `Allwmake`.
With a version that lacks some of the newer API, `Allwmake` builds a copy of the
code adjusted to it, `build/OpenFOAM-<version>`; `bin/adjustCodeToLoadedOFVersion -list`
lists the adjustments. The copy links to the files of the code, except a few
adjusted copies, listed in its `README`, which `Allwmake` remakes on each run;
`wmake` in a directory of the code builds it unadjusted, which fails with the older
version. Before v2506, interface positions are not transformed across
cyclic patches in plicRDF, so plicRDF stops on meshes with cyclic patches; the
periodic disc (`run/periodicDisc`) then needs v2506 or later.

The cases are set up and changed with
[dictator](https://github.com/isoAdvector/dictator), which rewrites only the value it
sets and keeps the comments, layout and full precision of the dictionary. Source it
(`source dictator/dictator`, for example in `~/.bashrc`) before running the cases;
the scripts inherit it from the shell they are started from. Two tools are
installed separately, each with its own `Allwmake`:

- [isoAdvector/streamFunctionWave](https://github.com/isoAdvector/streamFunctionWave):
  `streamFunFoam` and `calcStreamFunctionWave`, which make the initial and exact
  wave of the progressive wave cases (`run/progressiveWave`).
- [isoAdvector/surfacePosition](https://github.com/isoAdvector/surfacePosition):
  the function object `surfacePosition`, the wave gauge (`surfacePosition1`) of the
  progressive and standing wave cases (`run/progressiveWave`, `run/standingWave`).
  It writes the position of the interface at given locations, measured along a given
  direction (here that of gravity), every time step. A case with the gauge stops at the start if the
  library is not installed (`errors strict`).

The scans use `parmScanner` from
[isoAdvector/cfdTools](https://github.com/isoAdvector/cfdTools). The figure
scripts need python3 with numpy, scipy and matplotlib, and ParaView (pvbatch)
for the ParaView figures.

## Reproducing the papers

The runs of both papers were made with OpenFOAM v2512. Every figure, table and number of
the papers comes from a script in this repository. The summaries the scripts read
(`run/*/data/`) and the numbers of every table (`papers/part*/tab_*/*.txt`) are tracked
in git, so after a rerun `git status` and `git diff` show what changed. The full run data
take about 150 GB.

1. Install OpenFOAM v2512, this repository (`./Allwmake`), streamFunctionWave and
   surfacePosition (each with its own `./Allwmake`), dictator and cfdTools (Building),
   python3 with numpy, scipy and matplotlib, and ParaView (pvbatch).
   `run/quickTest/Allrun` checks the installation in about a minute.

2. Part I. Run the studies, each from `run/`:

   ```shell
   tiltedBox/Allrun
   skewedBox/Allrun
   refinedBox/Allrun
   standingWave/Allrun
   progressiveWave/Allrun        # several days on 32 cores
   ```

   `progressiveWave/Allrun` also runs the studies with `ddtCorr` off (Appendix A). The
   wall-clock times of the progressive wave come from separate timing runs on an
   otherwise idle machine: `python3 timing.py make`, `run` and `analyse` in
   `progressiveWave/timing_partI`, whose README gives the conditions.

3. Part II. Its fluxinator cases are made from the progressive-wave cases of Part I, so
   `progressiveWave/Allrun` must have run first. Then, from `run/`:

   ```shell
   taylorGreenVortex/Allrun
   prescribedFlows/Allrun
   periodicDisc/Allrun
   progressiveWave/AllrunPartII  # several days; ends with the timing runs
   ```

   The timing runs (`progressiveWave/timing_maxCo0.1_noGauge`) need cores with nothing
   else running; its README gives the conditions.

4. Make the figures and tables, in `papers/part1` for Part I and in `papers/part2` for
   Part II (there first `python3 browserData.py`, which the selection figures read):

   ```shell
   for d in fig_* tab_*
   do
       (
           cd "$d" || exit
           [ -f snapshots.py ] && pvbatch --force-offscreen-rendering snapshots.py
           if [ -f make.sh ]
           then
               bash make.sh
           elif grep -q paraview make.py
           then
               pvbatch --force-offscreen-rendering make.py
           else
               python3 make.py
           fi
       )
   done
   ```

   Each directory makes the figure or table with the same LaTeX label, whose caption in
   the paper names the directory (`papers/README.md`).

## Running

`run/quickTest/Allrun` builds nothing and takes about a minute: it runs
fluxinator on the accelerated disc and on one period of the wave, and prints the results
next to the expected values. Each case in `run/` has an `Allrun` for its studies
(`run/README.md`); `MAXJOBS` limits the number of simultaneous runs. All scripts use
the OpenFOAM environment and dictator of the shell they are started from. A case made by a
script in `run/` lists the entries the script changed in `log.dictator` (a parmScanner case
in `parmScannerParms`).

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

\* fluxinator is not endorsed by Dr. Heinz Doofenshmirtz or affiliated with
Doofenshmirtz Evil Incorporated.
