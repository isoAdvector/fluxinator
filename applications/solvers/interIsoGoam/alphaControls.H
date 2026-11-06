// Modified code Copyright (C) 2026 Johan Roenby
// SPDX-License-Identifier: GPL-3.0-or-later

const dictionary& alphaControls = mesh.solverDict(alpha1.name());

label nAlphaSubCycles(alphaControls.get<label>("nAlphaSubCycles"));

// Weight given to the latest outer iteration value of U and phi (vs.
// their oldTime value) when blending the velocity/flux used to advect
// alpha1 on outer (PIMPLE) iterations after the first. A value of 1
// disables blending (pure latest value); 0 always advects using the
// start-of-time-step value. See alphaEqnSubCycle.H.
scalar oCorrBlend
(
    alphaControls.getOrDefault<scalar>("outerCorrBlending", 0.5)
);
