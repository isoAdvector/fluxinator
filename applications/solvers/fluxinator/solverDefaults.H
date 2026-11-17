// Copyright (C) 2026 Johan Roenby
// SPDX-License-Identifier: GPL-3.0-or-later

// Run-time defaults that differ between the two solvers built from this source.
// fluxinator has the recommended ones. interIsoGoam (../interIsoGoam) keeps
// the defaults of interIsoFoam, so that the stored cases of Parts I and II give
// the same results. An entry in the case always overrides them.

#ifdef INTERISOGOAM
    const scalar defaultConvectionNewTimeWeight = 1;
    const bool defaultDdtCorr = true;
    const char* const defaultGravityModel = "basic";
    const label defaultNInitialReconstructions = 0;
#else
    const scalar defaultConvectionNewTimeWeight = 0;
    const bool defaultDdtCorr = false;
    const char* const defaultGravityModel = "geometric";
    const label defaultNInitialReconstructions = 3;
#endif
