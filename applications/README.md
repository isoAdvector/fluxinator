# Applications

| Directory | Content |
|---|---|
| `solvers/interIsoGoam` | the solver; its header and `timeWeights.H` list the switches |
| `test/convectionTest` | prescribed-flow test of the momentum convection (Part II) |
| `test/ddtCorrSwitch` | stock interFoam and interIsoFoam with a PIMPLE `ddtCorr` switch (Part I, Appendix A) |
| `utilities/setHydrostaticPressure` | sets the hydrostatic initial `p_rgh` |
| `utilities/surfacePosition` | surface position utility (submodule) |

Build all of them with `../Allwmake` in an OpenFOAM v2512 environment.
