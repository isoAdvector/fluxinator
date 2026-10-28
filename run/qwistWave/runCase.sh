#!/bin/bash
# runCase.sh <outdir> <nx> <nz> <dt> <endTime> [key=value ...]
# Stream-function wave of Qwist & Christensen (2023, Coastal Eng. 179, 104227):
# H = 0.3 m, T = 2 s, d = 1 m, zero Eulerian mean current, L = 5.4086 m, in a
# periodic one-wavelength domain, 2 m high, fixed time step (baseCase).
# Keys: app (default interIsoGoam), gravityModel (geometric | basic; basic =
# interIsoFoam bitwise), ddt (Euler | backward | ...), PIMPLE.<key>=<value>,
# libs (extra libraries, comma separated), conv (div(rhoPhi,U) scheme; '#' = blank).
# After the run only U and alpha.water are kept in the time directories
# (the last time is kept complete).
set +u; source /home/roenby/OpenFOAM/OpenFOAM-v2512/etc/bashrc > /dev/null 2>&1
out=$1; nx=$2; nz=$3; dt=$4; tEnd=$5; shift 5
app=interIsoGoam; gravityModel=""; ddt=""; libs=""; conv=""
for kv in "$@"; do k=${kv%%=*}; v=${kv#*=}; case $k in *.*) ;; *) eval "$k=\"\${v//#/ }\"";; esac; done
BASE=$(readlink -f "$(dirname "$0")/baseCase")
rm -rf "$out"; mkdir -p "$out"; cp -r $BASE/{0.orig,constant,system,Allrun} "$out"/; cd "$out" || exit 1
sed -i "s/^nx .*/nx $nx;/; s/^nz .*/nz $nz;/" system/blockMeshDict
foamDictionary -entry application -set $app  system/controlDict.run > /dev/null 2>&1
foamDictionary -entry deltaT      -set $dt   system/controlDict.run > /dev/null 2>&1
foamDictionary -entry endTime     -set $tEnd system/controlDict.run > /dev/null 2>&1
[ -n "$gravityModel" ] && foamDictionary -entry gravityModel -set $gravityModel constant/gravityProperties > /dev/null 2>&1
[ -n "$libs" ] && foamDictionary -entry libs -set "( ${libs//,/ } )" system/controlDict.run > /dev/null 2>&1
[ -n "$conv" ] && foamDictionary -entry 'divSchemes.div(rhoPhi,U)' -set "Gauss $conv" system/fvSchemes > /dev/null 2>&1
[ -n "$ddt" ] && foamDictionary -entry ddtSchemes.default -set "$ddt" system/fvSchemes > /dev/null 2>&1
for kv in "$@"; do k=${kv%%=*}; v=${kv#*=}
  case $k in PIMPLE.*) foamDictionary -entry "$k" -set "$v" system/fvSolution > /dev/null;; esac; done
echo "$out: app $app nx $nx nz $nz dt $dt endTime $tEnd $*" > caseInfo
./Allrun > log.Allrun 2>&1
grep -q "^End" log.$app || echo "FAILED $out"
last=$(foamListTimes -latestTime 2>/dev/null | tail -1)
for d in $(foamListTimes 2>/dev/null); do
    [ "$d" = "$last" ] && continue
    find "$d" -maxdepth 1 -type f ! -name U ! -name alpha.water -delete
    rm -rf "$d/uniform"
done
