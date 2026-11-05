#!/bin/bash
# runCaseNoForce.sh <case> [conv=<div(rhoPhi,U) scheme, '#' = blank>] [PIMPLE.<key>=<value> ...]
[ -n "${WM_PROJECT_DIR:-}" ] || { echo "Source the OpenFOAM v2512 environment first" >&2; exit 1; }
out=$1; shift
BASE=$(readlink -f "$(dirname "$0")/baseCaseNoForce")
rm -rf "$out"; mkdir -p "$out"; cp -r $BASE/{0.orig,constant,system,Allrun} "$out"/; cd "$out" || exit 1
for kv in "$@"; do k=${kv%%=*}; v=${kv#*=}; v=${v//#/ }
  case $k in
    conv) foamDictionary -entry 'divSchemes.div(rhoPhi,U)' -set "Gauss $v" system/fvSchemes > /dev/null;;
    PIMPLE.*) foamDictionary -entry "$k" -set "$v" system/fvSolution > /dev/null;;
  esac
done
echo "$out $*" > caseInfo
./Allrun > log.Allrun 2>&1
grep -q "^End" log.interIsoGoam || echo "FAILED $out"
