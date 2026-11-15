#!/bin/sh
# Copyright (C) 2026 Johan Roenby
# SPDX-License-Identifier: GPL-3.0-or-later
# Timing study at max Co 0.1 without the wave gauge (surfacePosition1 enabled false in
# controlDict.run; made when the gauge was interfaceHeight, which rebuilds a mesh search
# tree every step; the gauge is now surfacePosition, whose cost is negligible). 36 runs to
# t = 6 s, at most 8 at a time (./Allrun, unpinned), largest meshes first, on an otherwise
# idle machine. Needs the OpenFOAM v2606 environment.
cd "${0%/*}" || exit
[ -d final_95x20_r1 ] || ./makeCases.sh
ls -d *_380x80_r? *_190x40_r? *_95x20_r? | xargs -P 8 -I{} sh -c 'cd {} && ./Allrun > /dev/null 2>&1'
