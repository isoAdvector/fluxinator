#!/bin/bash
# Copyright (C) 2026 Johan Roenby
# SPDX-License-Identifier: GPL-3.0-or-later
# runCase.sh <name> <zw> <app> <gravityModel> <nLayers> <corr|uncorr> <zero|hydro> <nNonOrth>
# One untilted refined-box case in waterLevel_10steps/<name>_z<zw>, 10 steps of 0.01 s.
[ -n "$WM_PROJECT_DIR" ] || { echo "Source the OpenFOAM environment first" >&2; exit 1; }
type dictator > /dev/null 2>&1 || { echo "Source dictator first (github.com/isoAdvector/dictator)" >&2; exit 1; }
H=$(dirname "$(readlink -f "$0")"); name=$1; zw=$2; app=$3; gm=$4; nL=$5; sch=$6; start=$7; nc=$8
c=$H/waterLevel_10steps/${name}_z$zw; rm -rf $c; mkdir -p $H/waterLevel_10steps; cp -r $H/baseCase $c; cd $c
set -e   # stop if dictator fails (a missing file or an ambiguous entry)
{
dictator system/setAlphaFieldDict origin "(0 0 $zw)"
dictator system/fvSchemes interfaceForces.gravity.model $gm
dictator system/fvSchemes interfaceForces.gravity.nLayers $nL
dictator system/controlDict application $app
dictator system/controlDict endTime 0.1
dictator system/controlDict writeInterval 1
dictator system/fvSolution PIMPLE.nNonOrthogonalCorrectors $nc
if [ $sch = uncorr ]; then dictator system/fvSchemes laplacianSchemes.default "Gauss linear uncorrected"
  dictator system/fvSchemes snGradSchemes.default uncorrected; fi
# interFoam does not reconstruct the interface, which the reconstructedInterface surface needs:
# disable the surfaces function objects
if [ $app = interFoam ]; then
  for e in $(listDictEntries system/controlDict); do
    f=${e#functions.}; f=${f%.type}
    [ "$e" = "functions.$f.type" ] && [ "$(dictator system/controlDict $e)" = surfaces ] && dictator system/controlDict functions.$f.enabled false
  done
fi
} > log.dictator
set +e
[ $start = hydro ] && export HYDRO=1
./Allrun > /dev/null 2>&1; grep -q "^End" log.$app || echo "FAILED $c"
