# Standing wave on a locally refined mesh: experiments

Short (1 to 20 time step) runs used to find where the remaining spurious
currents of the geometric gravity model come from on this case. The wave starts
at rest at maximum elevation, so the exact velocity after one step of 0.01 s is
only 0.0154 m/s; anything above that is error.

`runCase.sh <outdir> [key=value ...]` copies a prepared case (`$SWBASE`,
default `/tmp/swbase`), applies the settings and runs it. Keys: `model`,
`nLayers`, `snGrad`, `lap`, `nOrthCorr`, `nCorr`, `amplitude`, `period`,
`alphaOrigin`, `recon`, `ptol`, `prel`, `dt`, `endTime`, `wint`, `app`.

Prepare a base case with the mesh already built, e.g.

```sh
mkdir /tmp/swbase && cp -r ../baseCase/{0.orig,constant,system} /tmp/swbase/
cd /tmp/swbase && rm -rf constant/polyMesh 0
blockMesh && topoSet -dict system/topoSetDict1
refineMesh -dict system/refineMeshDict1 -overwrite
```

Because the error depends strongly on where the interface falls relative to the
cell faces, compare medians over an ensemble of interface offsets spanning one
coarse cell rather than single runs (see FINDINGS.md).

`foamio.py` reads the ascii fields; split the velocity by phase with
`alpha.water` to separate the water solution from the air error.

Caveat: `foamDictionary -set` expands `$nx`-style variables when it rewrites a
dictionary, so edit `blockMeshDict` block counts with `-entry blocks` rather
than `-entry nx`.
