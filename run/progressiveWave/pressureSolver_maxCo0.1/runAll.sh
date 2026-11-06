#!/bin/sh
# Runs all cases at once with ./Allrun (24 single-core runs). Needs OpenFOAM v2606.
cd "${0%/*}" || exit
for c in *_B[0-3]; do (cd $c && ./Allrun > /dev/null 2>&1) & done
wait
