#!/bin/bash
# Copyright (C) 2026 Johan Roenby
# SPDX-License-Identifier: GPL-3.0-or-later
# runCase.sh <outdir> [key=value ...]  (copy of runCase.sh)
# With PIMPLE ddtCorr ${DDTCORR:-no}; interFoam and interIsoFoam run as interFoamDdtCorr and
# interIsoFoamDdtCorr (applications/ddtCorrSwitch). Base cases: baseUniformMesh, baseRefinedMesh.
# Keys: model nLayers snGrad laplacian nOrthCorr nCorr amplitude endTime app alphaOrigin
set +u
[ -n "${WM_PROJECT_DIR:-}" ] || { echo "Source the OpenFOAM v2606 environment first" >&2; exit 1; }
BASE=${SWBASE:-$(cd "${0%/*}" && pwd)/baseRefinedMesh}
out=$1; shift
rm -rf "$out"; cp -r "$BASE" "$out"
cd "$out" || exit 1
model=geometric; nLayers=2; snGrad=corrected; lap="Gauss linear corrected"
nOrthCorr=2; nCorr=3; amplitude=0.05; endTime=0.2; app=interIsoGoam; alphaOrigin=""; recon=isoAlpha; ptol=1e-07; prel=0.05; wint=20; dt=0.01; period=2
for kv in "$@"; do k=${kv%%=*}; v=${kv#*=}; eval "$k=\"\$v\""; done
foamDictionary -entry gravityModel -set "$model"   constant/gravityProperties > /dev/null
foamDictionary -entry nLayers      -set "$nLayers" constant/gravityProperties > /dev/null
foamDictionary -entry 'snGradSchemes.default'    -set "$snGrad" system/fvSchemes > /dev/null
foamDictionary -entry 'laplacianSchemes.default' -set "$lap"    system/fvSchemes > /dev/null
foamDictionary -entry 'PIMPLE.nNonOrthogonalCorrectors' -set "$nOrthCorr" system/fvSolution > /dev/null
foamDictionary -entry 'PIMPLE.nCorrectors'              -set "$nCorr"     system/fvSolution > /dev/null
foamDictionary -entry amplitude -set "$amplitude" system/setAlphaFieldDict > /dev/null
foamDictionary -entry period -set "$period" system/setAlphaFieldDict > /dev/null
[ -n "$alphaOrigin" ] && foamDictionary -entry origin -set "$alphaOrigin" system/setAlphaFieldDict > /dev/null
foamDictionary -entry 'solvers."alpha.*".reconstructionScheme' -set "$recon" system/fvSolution > /dev/null
foamDictionary -entry 'solvers.p_rgh.tolerance' -set "$ptol" system/fvSolution > /dev/null
foamDictionary -entry 'solvers.p_rgh.relTol'    -set "$prel" system/fvSolution > /dev/null
foamDictionary -entry deltaT -set "$dt" system/controlDict > /dev/null
foamDictionary -entry writeInterval -set "$wint" system/controlDict > /dev/null
foamDictionary -entry endTime   -set "$endTime"   system/controlDict > /dev/null
run=$app; case $app in interFoam|interIsoFoam) run=${app}DdtCorr;; esac
foamDictionary -entry application -set "$run"     system/controlDict > /dev/null
foamDictionary -entry PIMPLE.ddtCorr -set ${DDTCORR:-no} system/fvSolution > /dev/null
rm -rf 0 && cp -r 0.orig 0
setAlphaField > log.setAlphaField 2>&1
$run > log.run 2>&1
