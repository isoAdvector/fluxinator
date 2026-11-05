# Disc in a periodic box, full solver (Part II, Section 4)

The OpenFOAM tutorial `discInConstantFlowCyclicBCs`: a water disc of radius
0.15 m in air in a fully periodic box of 1 x 1 m with 50 x 50 cells, no gravity,
initial velocity U0 = (0.1, 0, 0.05) m/s, dt 0.005 s, ddtCorr off, run with
interIsoGoam.

| Cases | Content |
|---|---|
| `<convection>_rho{New,Old}` | uniform body force rho a, a = (0.5, 0, 0) m/s^2, with the density at the new or the old time level (`system/fvOptions`); exact U = U0 + a t, uniform p_rgh. `runCase.sh`, `baseCase` |
| `noForce_<convection>` | no force: the uniform velocity must be preserved; density ratio 1e3, 1 s. `runCaseNoForce.sh`, `baseCaseNoForce` |
| `noForceRatio1e4_<convection>` | the same with density ratio 1e4, 3 s. `runCaseNoForceRatio1e4.sh`, `baseCaseNoForceRatio1e4` |

`<convection>`: implicit or explicit convection with limitedLinearV (no suffix),
interfaceUpwind (IU), laxWendroff (LW), both (LWIU), upwind, linear, QUICK or
cubic face values. `writeUerr.py` writes the spurious velocity U - U0 - a t into
the time directories for viewing. The paper figure and table are made in
`papers/part2/fig_fullSolver` and `tab_fullSolver`.
