#!/bin/bash
# Pressure-solver settings study at max Co 0.1 (progressive wave, 30 s).
# Base set-ups: ../timing_maxCo0.1_t6/<solver>_<mesh>_r1 (themselves copies of
# solvers_maxCo_ddtCorrOff and finalSolver_maxCo) with endTime 30. Variants:
#   B0  as in the papers: PCG/DIC, nNonOrthogonalCorrectors 1, p_rghFinal 1e-12
#   B1  nNonOrthogonalCorrectors 0
#   B2  B1 + GAMG (smoother DIC) for p_rgh and p_rghFinal
#   B3  B2 + p_rghFinal tolerance 1e-9
# Needs dictator (setDictEntry/deleteDictEntry).
cd "${0%/*}" || exit
source ~/Projects/interFlow/dictator/dictator
for s in interFoam interIsoGoam final; do
  for m in 190x40 380x80; do
    for v in B0 B1 B2 B3; do
      c=${s}_${m}_$v; f=$c/system/fvSolution
      mkdir -p $c
      rsync -a --exclude polyMesh ../timing_maxCo0.1_t6/${s}_${m}_r1/{Allrun,0.orig,constant,system} $c/
      setDictEntry $c/system/controlDict.run endTime 30
      [ $v != B0 ] && setDictEntry $f PIMPLE.nNonOrthogonalCorrectors 0
      if [ $v = B2 ] || [ $v = B3 ]; then
        for p in p_rgh p_rghFinal; do
          setDictEntry $f solvers.$p.solver GAMG
          deleteDictEntry $f solvers.$p.preconditioner
          setDictEntry -add $f solvers.$p.smoother DIC
        done
      fi
      [ $v = B3 ] && setDictEntry $f solvers.p_rghFinal.tolerance 1e-09
    done
  done
done
