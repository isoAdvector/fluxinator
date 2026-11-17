# Applications

| Directory | Content |
|---|---|
| `solvers/fluxinator` | the solver to use, with the recommended defaults (`solverDefaults.H`); its header and `readControls.H` list the settings |
| `solvers/interIsoGoam` | the solver of the papers: the source of `fluxinator` with the defaults of interIsoFoam |
| `test/convectionTest` | prescribed-flow test of the momentum convection (Part II) |
| `test/ddtCorrSwitch` | stock interFoam and interIsoFoam with a PIMPLE `ddtCorr` switch (Part I, Appendix A) |
| `utilities/setHydrostaticPressure` | sets the hydrostatic initial `p_rgh` |

Build all of them with `../Allwmake` in an OpenFOAM v2606 environment.
