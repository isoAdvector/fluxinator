#!/bin/bash
# Runs a study again with PIMPLE ddtCorr no, as <study>_ddtCorrOff with the same case numbers
# (runStoredDdtCorrOff.sh). The study's cases must exist; their results are not needed. Finished
# cases are skipped. The 760x160 wave cases run on 4 subdomains, all others serially.
# Usage: runDdtCorrOffStudy.sh <study> [case ...]   (default: all cases; MAXJOBS: default nproc)
[ -n "${WM_PROJECT_DIR:-}" ] || { echo "Source the OpenFOAM v2512 environment first" >&2; exit 1; }
script=$(dirname "$(readlink -f "$0")")/runStoredDdtCorrOff.sh
study=$1; shift
[ $# -gt 0 ] || set -- $(ls "$study" | grep -E '^[0-9]{4}$')
mkdir -p "${study}_ddtCorrOff"
for c in "$@"; do
    s=$study/$c; o=${study}_ddtCorrOff/$c
    grep -qE "application: *streamFunFoam" "$s/parmScannerParms" 2>/dev/null && continue
    [ -e "$o/ddtCorrOff.info" ] && continue
    np=1; grep -q "nx: 760" "$s/parmScannerParms" 2>/dev/null && np=4
    printf '%s %s %s\n' "$s" "$o" "$np"
done | xargs -r -P "${MAXJOBS:-$(nproc)}" -L 1 "$script"
