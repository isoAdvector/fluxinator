# Shared modules and helpers

| File | Use |
|---|---|
| `kinematics.py` | stream function wave: mesh, exact wave, crest-column snapshots and errors |
| `hmatch.py` | wave height and celerity against the input wave (`analyse_input`) |
| `foamio.py` | reads ascii OpenFOAM fields |
| `summariseScan.py` | summary of a parmScanner study of the cases at rest (writes `data/<study>.txt`) |
| `snGradRhoReach` | function-object dictionary: the cells reached by a non-zero snGrad(rho) |
| `runStoredDdtCorrOff.sh` | reruns a stored case with PIMPLE `ddtCorr no` (the `_ddtCorrOff` studies) |
