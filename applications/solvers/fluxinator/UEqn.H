// Modified code Copyright (C) 2026 Johan Roenby
// SPDX-License-Identifier: GPL-3.0-or-later

    MRF.correctBoundaryVelocity(U);

    // Time derivative and convection w*div(rhoPhi, U^{n+1})
    // + (1 - w)*div(rhoPhi, U^n), with w = convectionNewTimeWeight
    // (readControls.H) and U^n = U.oldTime() also on later outer iterations.
    // Explicit convection (w 0) uses the mass flux rhoPhi (from phi^n) and the
    // convected velocity of the same time level; the implicit operator, zero
    // then, is not assembled.
    tmp<fvVectorMatrix> tddtConvection
    (
        convectionWeight == 0
      ? fvm::ddt(rho, U) + fvc::div(rhoPhi, U.oldTime(), "div(rhoPhi,U)")
      : convectionWeight < 1
      ? fvm::ddt(rho, U)
      + convectionWeight*fvm::div(rhoPhi, U)
      + (1 - convectionWeight)*fvc::div(rhoPhi, U.oldTime(), "div(rhoPhi,U)")
      : fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
    );

    fvVectorMatrix UEqn
    (
        tddtConvection
      + MRF.DDt(rho, U)
      + turbulence->divDevRhoReff(rho, U)
     ==
        fvOptions(rho, U)
    );

    #include "UEqnAddPorosity.H"

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                    mixture.surfaceTensionForce()
                  + gravityF
                  - fvc::snGrad(p_rgh)
                ) * mesh.magSf()
            )
        );

        fvOptions.correct(U);
    }
