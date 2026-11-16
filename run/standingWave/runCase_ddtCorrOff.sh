#!/bin/bash
# Copyright (C) 2026 Johan Roenby
# SPDX-License-Identifier: GPL-3.0-or-later
# runCase.sh <outdir> [key=value ...]  (copy of runCase.sh)
# With PIMPLE ddtCorr ${DDTCORR:-no}; interFoam and interIsoFoam run as interFoamDdtCorr and
# interIsoFoamDdtCorr (applications/ddtCorrSwitch). Base cases: baseUniformMesh, baseRefinedMesh.
# Keys: model nLayers snGrad lap nOrthCorr nCorr amplitude period alphaOrigin recon ptol prel dt wint endTime app
set +u
[ -n "${WM_PROJECT_DIR:-}" ] || { echo "Source the OpenFOAM v2606 environment first" >&2; exit 1; }
type dictator > /dev/null 2>&1 || { echo "Source dictator first (github.com/isoAdvector/dictator)" >&2; exit 1; }
BASE=${SWBASE:-$(cd "${0%/*}" && pwd)/baseRefinedMesh}
out=$1; shift
rm -rf "$out"; cp -r "$BASE" "$out"
cd "$out" || exit 1
model=geometric; nLayers=2; snGrad=corrected; lap="Gauss linear corrected"
nOrthCorr=2; nCorr=3; amplitude=0.05; endTime=0.2; app=interIsoGoam; alphaOrigin=""; recon=isoAlpha; ptol=1e-07; prel=0.05; wint=20; dt=0.01; period=2
for kv in "$@"; do k=${kv%%=*}; v=${kv#*=}; eval "$k=\"\$v\""; done
set -e   # stop if dictator fails (a missing file or an ambiguous entry)
{
dictator system/fvSchemes interfaceForces.gravity.model "$model"
dictator system/fvSchemes interfaceForces.gravity.nLayers "$nLayers"
dictator system/fvSchemes 'snGradSchemes.default' "$snGrad"
dictator system/fvSchemes 'laplacianSchemes.default' "$lap"
dictator system/fvSolution 'PIMPLE.nNonOrthogonalCorrectors' "$nOrthCorr"
dictator system/fvSolution 'PIMPLE.nCorrectors' "$nCorr"
dictator system/setAlphaFieldDict amplitude "$amplitude"
dictator system/setAlphaFieldDict period "$period"
[ -n "$alphaOrigin" ] && dictator system/setAlphaFieldDict origin "$alphaOrigin"
dictator system/fvSolution 'solvers."alpha.*".reconstructionScheme' "$recon"
dictator system/fvSolution 'solvers.p_rgh.tolerance' "$ptol"
dictator system/fvSolution 'solvers.p_rgh.relTol' "$prel"
dictator system/controlDict deltaT "$dt"
dictator system/controlDict writeInterval "$wint"
dictator system/controlDict endTime "$endTime"
run=$app; case $app in interFoam|interIsoFoam) run=${app}DdtCorr;; esac
dictator system/controlDict application "$run"
dictator system/fvSolution PIMPLE.ddtCorr ${DDTCORR:-no}
} > log.dictator
set +e
rm -rf 0 && cp -r 0.orig 0
setAlphaField > log.setAlphaField 2>&1
$run > log.run 2>&1
