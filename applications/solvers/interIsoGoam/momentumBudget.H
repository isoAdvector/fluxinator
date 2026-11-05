// Diagnostic only (no effect on the solution): split the velocity change of
// the last time step into accelerations, written at write times.
//
// PIMPLE entry:  momentumBudget  false;   (default off)
//
// With Euler and the isoAdvector mass flux, (rho - rho0)/dt = -div(rhoPhi), so
//   rho (U - U0)/dt = U0 div(rhoPhi) - div(rhoPhi, U) + F_reconstructed
// Fields [m/s^2]:
//   budget_dUdt  (U - U0)/dt                                   (total)
//   budget_conv  (U0 div(rhoPhi) - div(rhoPhi, U))/rho        (advective form,
//                same scheme as the momentum equation, final U)
//   budget_grav  reconstruct(gravity force |S|)/rho
//   budget_pres  -reconstruct(snGrad(p_rgh) |S|)/rho
//   budget_rest  dUdt - conv - grav - pres  (implicit-convection diagonal in
//                rAU, corrector lag, ...)
// rho = rho^{n+1}.
if
(
    pimple.dict().getOrDefault<bool>("momentumBudget", false)
 && runTime.writeTime()
)
{
    const dimensionedScalar dt(runTime.deltaT());

    volVectorField bdUdt("budget_dUdt", (U - U.oldTime())/dt);
    volVectorField bConv
    (
        "budget_conv",
        (U.oldTime()*fvc::div(rhoPhi) - fvc::div(rhoPhi, U, "div(rhoPhi,U)"))
       /rho
    );
    volVectorField bGrav
    (
        "budget_grav",
        fvc::reconstruct(gravity->force()*mesh.magSf())/rho
    );
    volVectorField bPres
    (
        "budget_pres",
        -fvc::reconstruct(fvc::snGrad(p_rgh)*mesh.magSf())/rho
    );
    volVectorField bRest("budget_rest", bdUdt - bConv - bGrav - bPres);

    bdUdt.write();
    bConv.write();
    bGrav.write();
    bPres.write();
    bRest.write();

    // Face ingredients of the convection term for a per-face split:
    // face velocity of the momentum convection scheme, mass flux, and
    // isoAdvector's water volume flux of the step (air mass flux =
    // rhoPhi - rho1*alphaPhi)
    tmp<fv::convectionScheme<vector>> tcs
    (
        fv::convectionScheme<vector>::New
        (
            mesh, rhoPhi, mesh.divScheme("div(rhoPhi,U)")
        )
    );
    // Face velocity as used by the momentum equation: theta-weighted
    // between U^{n+1} and U^n (see UEqn.H)
    const scalar bTheta = convTheta;
    surfaceVectorField
    (
        "budget_Uf",
        bTheta*tcs().interpolate(rhoPhi, U)
      + (1 - bTheta)*tcs().interpolate(rhoPhi, U.oldTime())
    ).write();
    surfaceScalarField("budget_rhoPhi", rhoPhi).write();
    surfaceScalarField("budget_alphaPhi", advector.alphaPhi()).write();
}
