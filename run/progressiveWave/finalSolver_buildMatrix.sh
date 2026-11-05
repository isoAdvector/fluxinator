#!/bin/bash
# Builds the Part 2 final-solver test matrix from the Part 1 interIsoGoam cases
# (solvers_fixedDt and solvers_maxCo), keeping the case numbers.
# Changes relative to the Part 1 case: explicit convection,
# div(rhoPhi,U) Gauss laxWendroff U phi interfaceUpwind alpha.water 1e-6 limitedLinearV 1,
# ddtCorr no; fields written every 0.2 s; the 760x160 cases decomposed into 4 subdomains.
# Usage: ./finalSolver_buildMatrix.sh   (creates finalSolver_fixedDt/ and finalSolver_maxCo/)
set -e
[ -n "${WM_PROJECT_DIR:-}" ] || { echo "Source the OpenFOAM v2512 environment first" >&2; exit 1; }
here=$(pwd)
for pair in "solvers_fixedDt:finalSolver_fixedDt" "solvers_maxCo:finalSolver_maxCo"; do
    src=${pair%%:*}; dst=${pair##*:}
    mkdir -p $dst
    for c in $(ls $src | grep -E '^[0-9]{4}$'); do
        grep -q "application: interIsoGoam " $src/$c/parmScannerParms 2>/dev/null || \
            grep -q "application: interIsoGoam$" $src/$c/parmScannerParms 2>/dev/null || continue
        d=$dst/$c; rm -rf $d; mkdir -p $d/constant
        cp -r $src/$c/0.orig $src/$c/system $src/$c/Allrun $d/
        for f in $(ls $src/$c/constant | grep -v polyMesh); do cp -r $src/$c/constant/$f $d/constant/; done
        cp baseCase_fixedDt/AllrunParallel $d/
        sed "s|$src/$c|$dst/$c|" $src/$c/parmScannerParms > $d/parmScannerParms
        cd $d
        foamDictionary -entry libs -set "(interfaceSchemes)" system/controlDict.run > /dev/null
        foamDictionary -entry writeControl -set adjustable system/controlDict.run > /dev/null
        foamDictionary -entry writeInterval -set 0.2 system/controlDict.run > /dev/null
        foamDictionary -entry 'divSchemes.div(rhoPhi,U)' -set 'Gauss laxWendroff U phi interfaceUpwind alpha.water 1e-6 limitedLinearV 1' system/fvSchemes > /dev/null
        foamDictionary -entry PIMPLE.explicitConvection -set yes system/fvSolution > /dev/null
        foamDictionary -entry PIMPLE.ddtCorr -set no system/fvSolution > /dev/null
        if grep -q "nx: 760" parmScannerParms; then
            foamDictionary -entry numberOfSubdomains -set 4 system/decomposeParDict > /dev/null
            foamDictionary -entry method -set simple system/decomposeParDict > /dev/null
            foamDictionary -entry simpleCoeffs.n -set "(4 1 1)" system/decomposeParDict > /dev/null
        fi
        cd $here
    done
done
ls finalSolver_fixedDt finalSolver_maxCo | tr '\n' ' '
