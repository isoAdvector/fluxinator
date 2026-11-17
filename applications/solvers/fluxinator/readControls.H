// Copyright (C) 2026 Johan Roenby
// SPDX-License-Identifier: GPL-3.0-or-later

// Run-time controls, read at the start (createFields.H) and once per time step.
//
// Time weights, each the weight of the value at the end of the time step, and
// 1 minus it the weight of the value at its start:
//
//   convectionNewTimeWeight (PIMPLE), w_c: the convected velocity,
//       w_c*div(rhoPhi, U^{n+1}) + (1 - w_c)*div(rhoPhi, U^n) (UEqn.H)
//   gravityNewTimeWeight (PIMPLE), w_g: the gravity force,
//       w_g*F^{n+1} + (1 - w_g)*F^n (gravityForce.H)
//   outerNewTimeWeight (solver dictionary of alpha), w_o: on the outer
//       correctors after the first, the velocity and flux that advect alpha,
//       w_o*U^{n+1} + (1 - w_o)*U^n, with U^{n+1} the latest estimate
//       (alphaEqnSubCycle.H; the first outer corrector advects with U^n)
//
// The defaults depend on nOuterCorrectors (convection also on the solver, see
// solverDefaults.H):
//
//   nOuterCorrectors 1   w_c 0 (fluxinator) or 1 (interIsoGoam), w_g 1
//   nOuterCorrectors > 1 w_c 1/2, w_g 1/2, w_o 1/2
//
// w_c 1 and w_g 1 is interIsoFoam, and w_c 1/2 Crank-Nicolson convection. With
// one outer corrector U^{n+1} and F^{n+1} are not known before the pressure
// solution, so w_c < 1 makes the convection term (partly) explicit, and w_c 0
// gives the mass flux and the convected velocity the same time level (Part II).
// The defaults for nOuterCorrectors > 1 make the step time-centred: alpha is
// advected with (U^n + U^{n+1})/2, and convection and gravity are averages of
// the start and the end of the step. The three weights belong together: in
// tests, w_o 1/2 with w_c 0 and w_g 1 damped a progressive wave strongly.
// w_c is also set in the registered convectionNewTimeWeight, for the
// laxWendroff scheme and the convectiveMomentumFlux function object.
//
// PIMPLE ddtCorr: the old-time flux correction of the predicted face flux
// (pEqn.H), default in solverDefaults.H.

{
    const bool outerCorrectors = pimple.nCorrPIMPLE() > 1;

    const scalar convectionWeightPrevious = convectionWeight;
    const scalar gravityWeightPrevious = gravityWeight;
    const scalar outerWeightPrevious = outerWeight;
    const bool ddtCorrPrevious = ddtCorr;

    convectionWeight =
        pimple.dict().getOrDefault<scalar>
        (
            "convectionNewTimeWeight",
            outerCorrectors ? 0.5 : defaultConvectionNewTimeWeight
        );

    gravityWeight =
        pimple.dict().getOrDefault<scalar>
        (
            "gravityNewTimeWeight",
            outerCorrectors ? 0.5 : 1
        );

    outerWeight =
        mesh.solverDict(alpha1.name()).getOrDefault<scalar>
        (
            "outerNewTimeWeight",
            0.5
        );

    ddtCorr = pimple.dict().getOrDefault<bool>("ddtCorr", defaultDdtCorr);

    convectionNewTimeWeight.value() = convectionWeight;

    if
    (
        convectionWeight != convectionWeightPrevious
     || gravityWeight != gravityWeightPrevious
     || (outerCorrectors && outerWeight != outerWeightPrevious)
     || ddtCorr != ddtCorrPrevious
    )
    {
        Info<< "Run-time controls: convectionNewTimeWeight " << convectionWeight
            << ", gravityNewTimeWeight " << gravityWeight;
        if (outerCorrectors)
        {
            Info<< ", outerNewTimeWeight " << outerWeight;
        }
        Info<< ", ddtCorr " << (ddtCorr ? "yes" : "no") << nl << endl;
    }

    // Explicit convection amplifies velocity differences at the interface
    // unless the convected face velocity there is taken from the upwind cell
    if (convectionWeight == 0 && !convectionSchemeChecked)
    {
        convectionSchemeChecked = true;
        bool hasInterfaceUpwind = false;
        for (const token& tok : mesh.divScheme("div(rhoPhi,U)"))
        {
            hasInterfaceUpwind =
                hasInterfaceUpwind
             || (tok.isWord() && tok.wordToken() == "interfaceUpwind");
        }
        if (!hasInterfaceUpwind)
        {
            WarningInFunction
                << "Explicit convection without the interfaceUpwind scheme in"
                << " div(rhoPhi,U): velocity differences are amplified in"
                << " interface cells. Use, for example," << nl
                << "    div(rhoPhi,U) Gauss laxWendroff interfaceUpwind"
                << " limitedLinearV 1;" << nl << endl;
        }
    }
}
