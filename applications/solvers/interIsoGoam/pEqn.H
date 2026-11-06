// Modified code Copyright (C) 2026 Johan Roenby
// SPDX-License-Identifier: GPL-3.0-or-later

{
    if (correctPhi)
    {
        rAU.ref() = 1.0/UEqn.A();
    }
    else
    {
        rAU = 1.0/UEqn.A();
    }

    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU()));
    volVectorField HbyA(constrainHbyA(rAU()*UEqn.H(), U, p_rgh));
    surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));

    // Old-time flux correction (PIMPLE entry ddtCorr, default yes, in
    // fluxinator no, see solverDefaults.H); with no, phiHbyA is the flux of
    // HbyA alone
    if (pimple.dict().getOrDefault<bool>("ddtCorr", defaultDdtCorr))
    {
        phiHbyA +=
            MRF.zeroFilter(fvc::interpolate(rho*rAU())*fvc::ddtCorr(U, phi, Uf));
    }
    MRF.makeRelative(phiHbyA);

    if (p_rgh.needReference())
    {
        fvc::makeRelative(phiHbyA, U);
        adjustPhi(phiHbyA, U, p_rgh);
        fvc::makeAbsolute(phiHbyA, U);
    }

    surfaceScalarField gravityF("gravityF", gravity->force());

    // Time weighting of the gravity force (gravityTheta, timeWeights.H):
    //   F = theta*F(rho^{n+1}, G^{n+1}) + (1 - theta)*F(rho^n, G^n), with
    //   F(rho^n, G^n) the force of the last pressure solution of the previous
    //   time step, kept fixed over the outer correctors. theta = 1 is the
    //   current force. Only the pressure equation is weighted
    //   (momentumPredictor no).
    if (gravityTheta != 1)
    {
        static label gravityOldIndex = runTime.timeIndex();

        if (!mesh.foundObject<surfaceScalarField>("gravityFLatest"))
        {
            for (const word name : {"gravityFLatest", "gravityFOld"})
            {
                auto* ptr = new surfaceScalarField
                (
                    IOobject
                    (
                        name,
                        runTime.timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE,
                        IOobject::REGISTER
                    ),
                    gravityF
                );
                ptr->store();
            }
        }

        surfaceScalarField& FLatest =
            mesh.lookupObjectRef<surfaceScalarField>("gravityFLatest");
        surfaceScalarField& FOld =
            mesh.lookupObjectRef<surfaceScalarField>("gravityFOld");

        if (runTime.timeIndex() != gravityOldIndex)
        {
            FOld = FLatest;
            gravityOldIndex = runTime.timeIndex();
        }
        FLatest = gravityF;

        gravityF = gravityTheta*gravityF + (1.0 - gravityTheta)*FOld;
    }

    surfaceScalarField phig
    (
        (
            mixture.surfaceTensionForce()
          + gravityF
        )*rAUf*mesh.magSf()
    );

    phiHbyA += phig;

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p_rgh, U, phiHbyA, rAUf, MRF);

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqn
        (
            fvm::laplacian(rAUf, p_rgh) == fvc::div(phiHbyA)
        );

        p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));

        p_rghEqn.solve(p_rgh.select(pimple.finalInnerIter()));

        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA - p_rghEqn.flux();

            p_rgh.relax();

            U = HbyA + rAU()*fvc::reconstruct((phig - p_rghEqn.flux())/rAUf);
            U.correctBoundaryConditions();
            fvOptions.correct(U);
        }
    }

    #include "continuityErrs.H"

    // Correct Uf if the mesh is moving
    fvc::correctUf(Uf, U, phi);

    // Make the fluxes relative to the mesh motion
    fvc::makeRelative(phi, U);

    p == p_rgh + rho*gh;

    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
        p_rgh = p - rho*gh;
    }

    if (!correctPhi)
    {
        rAU.clear();
    }
}
