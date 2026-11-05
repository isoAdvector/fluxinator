# Libraries

| Directory | Library | Content |
|---|---|---|
| `geometricVoF` | `newGeometricVoF` | isoAdvector and interface reconstruction used by interIsoGoam |
| `gravityModels` | `gravityModels` | gravity models `basic`, `geometric` and `constantHeight` (Part I) |
| `interfaceSchemes` | `interfaceSchemes` | interpolation schemes `interfaceUpwind` and `laxWendroff` (Part II) |
| `functionObjects` | `interfaceFunctionObjects` | function object `convectiveMomentumFlux` (Part II) |

Build all of them with `../Allwmake` in an OpenFOAM v2512 environment. Example settings
are in the test cases in `../run`.
