// Run-time defaults that differ between the two solvers built from this source.
// interIsoGoam keeps the defaults of interIsoFoam, so that its stored cases give
// the same results; fluxinator (../fluxinator) has the recommended ones. An entry
// in the case always overrides them.

#ifdef FLUXINATOR
    const bool defaultExplicitConvection = true;
    const bool defaultDdtCorr = false;
    const char* const defaultGravityModel = "geometric";
#else
    const bool defaultExplicitConvection = false;
    const bool defaultDdtCorr = true;
    const char* const defaultGravityModel = "basic";
#endif
