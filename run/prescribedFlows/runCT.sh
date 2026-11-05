#!/bin/bash
# runCT.sh <case> <mode reset|evolve> <disc origin "x 0 z"> <endTime> "<flow dictionary entries>"
# REF=massCentroid sets the reference point (default cellCentre, see convectionTest.C)
# e.g. runCT.sh rot reset "0.5 0 0.75" 1 "Omega 2; rotationCentre (0.5 0 0.5);"
set +u; source /home/roenby/OpenFOAM/OpenFOAM-v2512/etc/bashrc > /dev/null 2>&1
out=$1; mode=$2; origin=$3; tEnd=$4; flow=$5
BASE=$(readlink -f "$(dirname "$0")/baseCase")
rm -rf "$out"; mkdir -p "$out"; cp -r $BASE/* "$out"/; cd "$out" || exit 1
[ -n "$N" ] && sed -i "s/^n 50;/n $N;/" system/blockMeshDict
[ -n "$DT" ] && sed -i "s/deltaT 0.005;/deltaT $DT;/" system/controlDict
sed -i "s/origin (0.5 0 0.75)/origin ($origin)/" system/setAlphaFieldDict
sed -i "s/endTime 1;/endTime $tEnd;/" system/controlDict
cat > system/convectionTestDict <<EOD
FoamFile { version 2.0; format ascii; class dictionary; object convectionTestDict; }
mode $mode;
${REF:+reference $REF;}
flow { $flow }
#include "variants"
EOD
./Allrun > log.Allrun 2>&1
grep -q "^End" log.convectionTest || echo "FAILED $out"
