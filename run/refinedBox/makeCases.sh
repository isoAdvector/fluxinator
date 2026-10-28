#!/bin/bash
# Prints the case list for runCase.sh: every configuration at every water level
# z_w = 0.48, 0.4825, ..., 0.52 (one coarse cell, 17 levels).
# Columns: name z_w application gravityModel nLayers schemes start nNonOrthogonalCorrectors
configs="
iF_corr  interFoam    basic          1 corr   zero  0
iF_nc2   interFoam    basic          1 corr   zero  2
iF_hyd   interFoam    basic          1 corr   hydro 0
iF_unc   interFoam    basic          1 uncorr zero  0
iIF_corr interIsoGoam basic          1 corr   zero  0
iIF_nc2  interIsoGoam basic          1 corr   zero  2
iIF_hyd  interIsoGoam basic          1 corr   hydro 0
iIF_unc  interIsoGoam basic          1 uncorr zero  0
cH_zero  interIsoGoam constantHeight 1 corr   zero  0
cH_nc2   interIsoGoam constantHeight 1 corr   zero  2
cH_hyd   interIsoGoam constantHeight 1 corr   hydro 0
cH_unc   interIsoGoam constantHeight 1 uncorr zero  0
g1_zero  interIsoGoam geometric      1 corr   zero  0
g1_nc2   interIsoGoam geometric      1 corr   zero  2
g1_hyd   interIsoGoam geometric      1 corr   hydro 0
g1_unc   interIsoGoam geometric      1 uncorr zero  0
g2_zero  interIsoGoam geometric      2 corr   zero  0
g2_nc2   interIsoGoam geometric      2 corr   zero  2
g2_hyd   interIsoGoam geometric      2 corr   hydro 0
g2_unc   interIsoGoam geometric      2 uncorr zero  0
g3_hyd   interIsoGoam geometric      3 corr   hydro 0"
for i in $(seq 0 16); do
    z=$(awk -v i=$i 'BEGIN{printf "%.4f", 0.48 + 0.0025*i}')
    echo "$configs" | awk -v z=$z 'NF {print $1, z, $2, $3, $4, $5, $6, $7}'
done
