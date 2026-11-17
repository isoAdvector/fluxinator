// Copyright (C) 2026 Johan Roenby
// SPDX-License-Identifier: GPL-3.0-or-later

// Face gravity force of this outer corrector, used by the momentum predictor
// and the pressure equation. Weighted in time with gravityNewTimeWeight w_g
// (readControls.H):
//
//   F = w_g*F(rho^{n+1}, G^{n+1}) + (1 - w_g)*F(rho^n, G^n)
//
// with F(rho^n, G^n) the force of the last outer corrector of the previous
// time step (the initial force in the first time step, or the current force
// when the weight changes from 1 during a run). w_g = 1 is the current force,
// and the forces are then not stored (createFields.H).

surfaceScalarField gravityF("gravityF", gravity->force());

if (gravityWeight != 1)
{
    if (!gravityFOldPtr)
    {
        // The weight was 1 until now: start from the current force
        gravityFLatestPtr.reset
        (
            new surfaceScalarField("gravityFLatest", gravityF)
        );
        gravityFOldPtr.reset(new surfaceScalarField("gravityFOld", gravityF));
    }
    else if (pimple.firstIter())
    {
        gravityFOldPtr() = gravityFLatestPtr();
    }
    gravityFLatestPtr() = gravityF;

    gravityF = gravityWeight*gravityF + (1.0 - gravityWeight)*gravityFOldPtr();
}
