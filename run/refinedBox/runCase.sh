#!/bin/bash
# Copyright (C) 2026 Johan Roenby
# SPDX-License-Identifier: GPL-3.0-or-later
# runCase.sh <name> <zw> <app> <gravityModel> <nLayers> <corr|uncorr> <zero|hydro> <nNonOrth>
# One untilted refined-box case in waterLevel_10steps/<name>_z<zw>, 10 steps of 0.01 s.
[ -n "$WM_PROJECT_DIR" ] || { echo "Source the OpenFOAM environment first" >&2; exit 1; }
H=$(dirname "$(readlink -f "$0")"); name=$1; zw=$2; app=$3; gm=$4; nL=$5; sch=$6; start=$7; nc=$8
c=$H/waterLevel_10steps/${name}_z$zw; rm -rf $c; mkdir -p $H/waterLevel_10steps; cp -r $H/baseCase $c; cd $c
foamDictionary -entry origin -set "(0 0 $zw)" system/setAlphaFieldDict > /dev/null 2>&1
foamDictionary -entry gravityModel -set $gm constant/gravityProperties > /dev/null 2>&1
foamDictionary -entry nLayers -set $nL constant/gravityProperties > /dev/null 2>&1
foamDictionary -entry application -set $app system/controlDict > /dev/null 2>&1
foamDictionary -entry endTime -set 0.1 system/controlDict > /dev/null 2>&1
foamDictionary -entry writeInterval -set 1 system/controlDict > /dev/null 2>&1
foamDictionary -entry PIMPLE.nNonOrthogonalCorrectors -set $nc system/fvSolution > /dev/null 2>&1
if [ $sch = uncorr ]; then foamDictionary -entry laplacianSchemes.default -set "Gauss linear uncorrected" system/fvSchemes > /dev/null 2>&1
  foamDictionary -entry snGradSchemes.default -set uncorrected system/fvSchemes > /dev/null 2>&1; fi
[ $app = interFoam ] && for f in $(foamDictionary -keywords -entry functions system/controlDict); do [ "$(foamDictionary -entry functions.$f.type -value system/controlDict)" = surfaces ] && foamDictionary -entry functions.$f -remove system/controlDict > /dev/null; done
[ $start = hydro ] && export HYDRO=1
./Allrun > /dev/null 2>&1; grep -q "^End" log.$app || echo "FAILED $c"
