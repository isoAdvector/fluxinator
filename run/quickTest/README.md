# Quick test

`./Allrun` in an OpenFOAM v2512 environment, after `../../Allwmake`, takes a few minutes on one
core. It runs the final solver of Part II on the accelerated disc (Section 4) and on one period of the
progressive wave (needs `streamFunFoam`, skipped otherwise), and prints the results next to the
expected values.

The solver can stop with a memory error at exit, after it has printed `End`. This comes from static
objects defined in two libraries and does not affect the results, which are written before it happens.
