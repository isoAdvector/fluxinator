#!/bin/sh
# Timing study at max Co 0.1 without the interfaceHeight gauge (enabled false in
# controlDict.run; it rebuilds a mesh search tree every step). 36 runs to t = 6 s,
# at most 8 at a time (./Allrun, unpinned), largest meshes first, on an otherwise idle
# machine. Needs the OpenFOAM v2606 environment.
cd "${0%/*}" || exit
[ -d final_95x20_r1 ] || ./makeCases.sh
ls -d *_380x80_r? *_190x40_r? *_95x20_r? | xargs -P 8 -I{} sh -c 'cd {} && ./Allrun > /dev/null 2>&1'
