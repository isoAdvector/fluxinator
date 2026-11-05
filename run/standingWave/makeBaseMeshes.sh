#!/bin/bash
# Builds the meshes of the prepared base cases copied by runCase.sh:
#   baseUniformMesh  blockMesh only
#   baseRefinedMesh  blockMesh, then one level of refinement around the interface
#                    (topoSetDict1, refineMeshDict1), as baseCase/Allrun
# Usage: ./makeBaseMeshes.sh  (OpenFOAM environment sourced)
cd "${0%/*}" || exit 1
[ -n "$WM_PROJECT_DIR" ] || { echo "Source the OpenFOAM environment first" >&2; exit 1; }
( cd baseUniformMesh && rm -rf constant/polyMesh && blockMesh > log.blockMesh 2>&1 )
( cd baseRefinedMesh && rm -rf constant/polyMesh && blockMesh > log.blockMesh 2>&1 \
  && topoSet -dict system/topoSetDict1 > log.topoSet 2>&1 \
  && refineMesh -dict system/refineMeshDict1 -overwrite > log.refineMesh 2>&1 )
echo "baseUniformMesh: $(grep -c . baseUniformMesh/constant/polyMesh/owner) lines in owner; baseRefinedMesh built"
