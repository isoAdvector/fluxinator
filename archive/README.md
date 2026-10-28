# Archive

Superseded code kept for the record. Nothing here is compiled.

- `gravityModels/geometric.v0`, `geometric.v1`, `geometric.v2`: earlier
  versions of the geometric gravity model. They were moved out of
  `gravityModels/` because their files share the name
  `geometricGravityModel.[CH]` with the current `gravityModels/geometric`, so
  `gravityModels/lnInclude` could link to any one of them.
- `gravityModels/gravityReconNoProjection`, `gravityReconOnlyInterfaceCells`:
  variants of `gravityRecon` that were never in `gravityModels/Make/files`.

Superseded test cases (e.g. `run/tiltedBoxRefined`, replaced by
`run/refinedBox`) were removed from the repository on 2026-10-28; they are in
the git history and, with their case data, in `../geometricGravity_archive/`
(see its `MOVES.txt`).
