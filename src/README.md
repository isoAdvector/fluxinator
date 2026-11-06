# Libraries

| Directory | Library | Content |
|---|---|---|
| `geometricVoF` | `newGeometricVoF` | isoAdvector and interface reconstruction used by interIsoGoam |
| `gravityModels` | `gravityModels` | gravity models `basic`, `geometric` and `constantHeight` (Part I) |
| `interfaceSchemes` | `interfaceSchemes` | interpolation schemes `interfaceUpwind` and `laxWendroff` (Part II) |
| `functionObjects` | `interfaceFunctionObjects` | function object `convectiveMomentumFlux` (Part II) |

`newGeometricVoF` is in the namespace `Foam::isoAdvection`, so that it can be loaded next to OpenFOAM's
`geometricVoF`, which many function objects pull in. Its sampled surface type is therefore
`reconstructedInterface`, not `interface`.

Build all of them with `../Allwmake` in an OpenFOAM v2512 environment. Example settings
are in the test cases in `../run`.
