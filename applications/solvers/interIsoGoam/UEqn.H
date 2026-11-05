    MRF.correctBoundaryVelocity(U);

    // Convection theta*div(rhoPhi, U^{n+1}) + (1 - theta)*div(rhoPhi, U^n),
    // U^n = U.oldTime() also on later outer iterations; convTheta and its
    // defaults (explicitConvection, nOuterCorrectors) in timeWeights.H.
    // Explicit convection (theta 0) uses the mass flux rhoPhi (from phi^n)
    // and the convected velocity of the same time level.

    fvVectorMatrix UEqn
    (
        convTheta < 1
      ? tmp<fvVectorMatrix>
        (
            fvm::ddt(rho, U)
          + convTheta*fvm::div(rhoPhi, U)
          + (1 - convTheta)*fvc::div(rhoPhi, U.oldTime(), "div(rhoPhi,U)")
          + MRF.DDt(rho, U)
          + turbulence->divDevRhoReff(rho, U)
         ==
            fvOptions(rho, U)
        )
      : tmp<fvVectorMatrix>
        (
            fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
          + MRF.DDt(rho, U)
          + turbulence->divDevRhoReff(rho, U)
         ==
            fvOptions(rho, U)
        )
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
                  + gravity->force()
                  - fvc::snGrad(p_rgh)
                ) * mesh.magSf()
            )
        );

        fvOptions.correct(U);
    }
