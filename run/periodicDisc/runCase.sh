#!/bin/bash
# Copyright (C) 2026 Johan Roenby
# SPDX-License-Identifier: GPL-3.0-or-later
# runCase.sh <case> [conv=<div(rhoPhi,U) scheme, '#' = blank>] [force=old] [PIMPLE.<key>=<value> ...]
# The body acceleration is set in baseCase/system/fvOptions (density new; force=old: old).
[ -n "${WM_PROJECT_DIR:-}" ] || { echo "Source the OpenFOAM v2606 environment first" >&2; exit 1; }
type dictator > /dev/null 2>&1 || { echo "Source dictator first (github.com/isoAdvector/dictator)" >&2; exit 1; }
out=$1; shift
BASE=$(readlink -f "$(dirname "$0")/baseCase")
rm -rf "$out"; mkdir -p "$out"; cp -r $BASE/{0.orig,constant,system,Allrun} "$out"/; cd "$out" || exit 1
set -e   # stop if dictator fails (a missing file or an ambiguous entry)
for kv in "$@"; do k=${kv%%=*}; v=${kv#*=}; v=${v//#/ }
  case $k in
    conv) dictator system/fvSchemes 'divSchemes.div(rhoPhi,U)' "Gauss $v";;
    force) [ "$v" = old ] && sed -i 's/eqn += rho\*/eqn += rho.oldTime()*/' system/fvOptions;;
    PIMPLE.*) dictator system/fvSolution "$k" "$v";;
  esac
done > log.dictator
set +e
./Allrun > log.Allrun 2>&1
grep -q "^End" log.interIsoGoam || echo "FAILED $out"
