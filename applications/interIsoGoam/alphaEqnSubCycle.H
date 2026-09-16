if (pimple.nCorrPIMPLE() > 1)
{
    // If nOuterCorrectors > 1 then for all but the first loop the advection
    // of alpha is done using a blend,
    // oCorrBlend*phiNew + (1-oCorrBlend)*phi,
    // where phi is the flux at the beginning of the time step and phiNew
    // is the flux estimate at the end of the time step from the previous
    // outer iteration. Similarly for U and UNew.
    if (pimple.firstIter())
    {
        // Storing the interface data (centre, normal, interfaceCell,
        // interfaceLabels) corresponding to alpha1.oldTime(), i.e. before
        // it gets advected below, so it can be cheaply restored on later
        // outer iterations instead of via a fresh, expensive
        // reconstruction.
        advector.surf().storeOld();
    }
    else
    {
        // Resetting alpha1 to its value at the start of the time step
        // before repeating the alpha advection in this outer iteration.
        alpha1 = alpha1.oldTime();

        // Restoring the interface data stored above so that it is
        // consistent with the alpha1 reset, without the cost of a fresh
        // reconstruction.
        advector.surf().resetToOld();

        // Storing the latest U and phi so they can be restored after the
        // alpha advection below, then temporarily overwriting them with
        // the blend used to advect the interface.
        U.storePrevIter();
        phi.storePrevIter();
        U = oCorrBlend*U + (1 - oCorrBlend)*U.oldTime();
        phi = oCorrBlend*phi + (1 - oCorrBlend)*phi.oldTime();
    }
}

if (nAlphaSubCycles > 1)
{
    dimensionedScalar totalDeltaT = runTime.deltaT();
    surfaceScalarField rhoPhiSum
    (
        mesh.newIOobject("rhoPhiSum"),
        mesh,
        dimensionedScalar(rhoPhi.dimensions(), Zero)
    );

    for
    (
        subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
        !(++alphaSubCycle).end();
    )
    {
        #include "alphaEqn.H"
        rhoPhiSum += (runTime.deltaT()/totalDeltaT)*rhoPhi;
    }

    rhoPhi = rhoPhiSum;
}
else
{
    #include "alphaEqn.H"
}

if (!pimple.firstIter())
{
    // Restoring U and phi to their value at the latest outer iteration,
    // stored above prior to blending.
    U = U.prevIter();
    phi = phi.prevIter();
}

rho == alpha1*rho1 + alpha2*rho2;
