#!/bin/bash
# Runs the given case directories, each with its own Allrun, at most MAXJOBS at a time
# (default: the number of cores). A case with a solver log has been run and is skipped.
# Usage: runCases.sh <solver log name> <case> ...
[ -n "${WM_PROJECT_DIR:-}" ] || { echo "Source the OpenFOAM v2512 environment first" >&2; exit 1; }
solverLog=$1; shift
for c in "$@"; do
    [ -e "$c/$solverLog" ] || printf '%s\n' "$c"
done | xargs -r -P "${MAXJOBS:-$(nproc)}" -I{} sh -c './{}/Allrun > {}/log.Allrun 2>&1 || echo "FAILED {}"'
for c in "$@"; do
    grep -q "^End" "$c/$solverLog" 2>/dev/null || echo "NOT FINISHED: $c"
done
