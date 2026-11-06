#!/bin/sh
# Runs the three cases at once with ./Allrun. Needs the OpenFOAM v2606 environment.
cd "${0%/*}" || exit
for c in final_190x40_tol*; do (cd $c && ./Allrun > /dev/null 2>&1) & done
wait
