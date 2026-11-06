#!/bin/bash
# runStored.sh <storedCase> <outdir> [nProcs]
# Rerun a stored Part I case (run/<case>/<study>/<run>, any Allrun-driven case)
# with the single change PIMPLE ddtCorr no (DDTCORR=yes in the environment for
# regression checks). Stock solvers are
# replaced by copies with the ddtCorr switch (applications/test/ddtCorrSwitch):
# interFoam -> interFoamDdtCorr, interIsoFoam -> interIsoFoamDdtCorr; with
# ddtCorr yes these reproduce the stored runs bitwise (reg/). Figure rendering
# (pvpython) is skipped. With nProcs > 1 the wave case runs with AllrunParallel
# on nProcs subdomains and is reconstructed afterwards.
[ -n "${WM_PROJECT_DIR:-}" ] || { echo "Source the OpenFOAM v2512 environment first" >&2; exit 1; }
S=$(readlink -f "$1"); out=$2; np=${3:-1}
rm -rf "$out"; mkdir -p "$out"
for f in $S/*; do
    b=$(basename $f)
    case $b in
        processor*|postProcessing|log.*|*.png|*.foam|dynamicCode) continue;;
    esac
    [[ $b =~ ^[0-9.e+-]+$ ]] && continue     # time directories
    cp -r $f "$out"/
done
rm -rf "$out"/constant/polyMesh
cd "$out" || exit 1
sed -i 's/^\( *\)pvpython/\1: pvpython/' Allrun
for d in system/controlDict system/controlDict.run; do
    [ -f $d ] || continue
    app=$(foamDictionary -entry application -value $d)
    case $app in
        interFoam)    foamDictionary -entry application -set interFoamDdtCorr $d > /dev/null;;
        interIsoFoam) foamDictionary -entry application -set interIsoFoamDdtCorr $d > /dev/null;;
    esac
done
foamDictionary -entry PIMPLE.ddtCorr -set ${DDTCORR:-no} system/fvSolution > /dev/null
echo "stored: $1  change: PIMPLE.ddtCorr ${DDTCORR:-no}  nProcs: $np" > ddtCorrOff.info
if [ "$np" -gt 1 ]; then
    cp $(dirname $(readlink -f "$0"))/../progressiveWave/baseCase_fixedDt/AllrunParallel .
    foamDictionary -entry numberOfSubdomains -set $np system/decomposeParDict > /dev/null
    foamDictionary -entry method -set simple system/decomposeParDict > /dev/null
    foamDictionary -entry simpleCoeffs.n -set "($np 1 1)" system/decomposeParDict > /dev/null
    ./AllrunParallel > log.Allrun 2>&1
    reconstructPar > log.reconstructPar 2>&1 && rm -rf processor*
else
    ./Allrun > log.Allrun 2>&1
fi
app=$(foamDictionary -entry application -value system/controlDict)
grep -q "^End" log.$app || echo "FAILED $out"
