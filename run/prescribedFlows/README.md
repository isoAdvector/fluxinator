# Prescribed flows (Part II, Sections 4 and 5)

A disc of water (radius 0.15 m, 1000 kg/m3) in air (1 kg/m3) in a 2-D box of
1 x 1 m with 50 x 50 cells. The velocity is prescribed (stream function), the
face fluxes are divergence free to round-off, and only the momentum equation is
solved, for several convection variants at once, by the test application
`applications/test/convectionTest` (see the header of `convectionTest.C`).

| Case | Flow | Mode |
|---|---|---|
| `uniform` | U0 = (0.2, 0, 0.05) m/s | one step |
| `acceleration`, `accelerationCarried` | U0 = (0.1, 0, 0.05) m/s + a t, a = (0.5, 0, 0) m/s^2 | one step, carried |
| `shear`, `shearCarried` | u = (z - 0.5, 0, -0.2), disc at (0.5, 0.75) | one step, carried |
| `rotation`, `rotationCarried` | solid body, Omega = 2 about (0.5, 0.5) | one step, carried |
| `<case>_n25`, `<case>_n100` | the same on 25 x 25 and 100 x 100 cells (dt scaled) | |

One step: every step starts from the exact cell-centre velocity (the error of
one step). Carried: the computed velocity is carried from step to step.

`runCT.sh <case> <mode oneStep|carried> "<disc origin>" <endTime> "<flow entries>"`
makes a case from `baseCase` and runs it; `summarise.sh <case>` prints the
largest errors per variant. The paper figures and tables are made in
`papers/part2` (`tab_rigUniform`, `tab_rigFaces`, `tab_rigConvergence`,
`fig_rigErrors`, `fig_rigCasesFlows`).
