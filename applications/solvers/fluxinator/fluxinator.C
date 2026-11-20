/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016 DHI
    Copyright (C) 2017 OpenCFD Ltd.
    Copyright (C) 2018, 2025-2026 Johan Roenby
    Copyright (C) 2019-2020 DLR
    Copyright (C) 2020 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    fluxinator

Group
    grpMultiphaseSolvers

Description
    Two-phase Volume-of-Fluid solver for two incompressible, isothermal
    immiscible fluids, derived from interFoam and interIsoFoam, with
    geometric interface advection (isoAdvector), a consistent gravity force
    (Part I) and a consistent momentum convection (Part II). Optional mesh
    motion and mesh topology changes including adaptive re-meshing.

    Defaults (solverDefaults.H), each overridden by an entry in the case:

    - interfaceForces.gravity.model geometric (system/fvSchemes)
    - PIMPLE convectionNewTimeWeight 0 (explicit convection)
    - PIMPLE ddtCorr no
    - nInitialReconstructions 3 (solver dictionary of alpha)

    The convection scheme is set in system/fvSchemes, as in any solver.
    Explicit convection requires the interfaceUpwind scheme, and the solver
    warns without it. For waves use

        div(rhoPhi,U) Gauss laxWendroff interfaceUpwind limitedLinearV 1;

    The same source builds interIsoGoam (../interIsoGoam), the solver of
    Parts I and II, which has the defaults of interIsoFoam (basic gravity
    model, convectionNewTimeWeight 1, ddtCorr yes, nInitialReconstructions
    0). With these defaults and one outer corrector it is interIsoFoam.

    Settings, all selected at run time:

    - Gravity model (interfaceForces.gravity.model in system/fvSchemes).
      The geometric model evaluates the gravity term at the height of the
      reconstructed interface (Part I). The basic model is that of interFoam
      and interIsoFoam.

    - Time weights (readControls.H), each the weight of the value at the end
      of the time step: convectionNewTimeWeight (PIMPLE) for the convected
      velocity of div(rhoPhi, U), 0 explicit, 1/2 Crank-Nicolson, 1
      implicit (Part II); gravityNewTimeWeight (PIMPLE) for the gravity
      force; outerNewTimeWeight (solver dictionary of alpha) for the
      velocity and flux that advect alpha on the outer correctors after the
      first. With nOuterCorrectors > 1 all three default to 1/2, which
      centres the step in time.

    - PIMPLE ddtCorr.

    - nInitialReconstructions (solver dictionary of alpha): reconstructions
      of the initial interface before the first time step.

    References:
    \verbatim
        Roenby, J., Missios, K., Møller, K., Aliyar, S. and Scheufler, H. (2026).
        Consistent interfacial flow numerics - Part I: The gravity force term.
        arXiv:2611.24841, https://arxiv.org/abs/2611.24841

        Roenby, J. (2026).
        Consistent interfacial flow numerics - Part II: The momentum
        convection term.
        arXiv:2611.24866, https://arxiv.org/abs/2611.24866

        Roenby, J., Bredmose, H. and Jasak, H. (2016).
        A computational method for sharp interface advection
        Royal Society Open Science, 3
        doi 10.1098/rsos.160405
    \endverbatim

    isoAdvector code supplied by Johan Roenby, STROMNING (2018)

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "isoAdvector.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "incompressibleInterPhaseTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "CorrectPhi.H"
#include "fvcSmooth.H"
#include "dynamicRefineFvMesh.H"
#include "gravityModel.H"
#include "uniformDimensionedFields.H"
#include "solverDefaults.H"

using namespace Foam::fluxinatorVoF;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Two-phase Volume-of-Fluid solver for two incompressible, isothermal"
        " immiscible fluids, with geometric interface advection (isoAdvector),"
        " a consistent gravity force and a consistent momentum convection."
        " Derived from interFoam and interIsoFoam."
    #ifdef INTERISOGOAM
        "\ninterIsoGoam is fluxinator with the defaults of interIsoFoam."
    #endif
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "initCorrectPhi.H"
    #include "createUfIfPresent.H"

    #include "porousCourantNo.H"
    #include "setInitialDeltaT.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readDyMControls.H"
        #include "porousCourantNo.H"
        #include "porousAlphaCourantNo.H"
        #include "setDeltaT.H"

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "readControls.H"

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            if (pimple.firstIter() || moveMeshOuterCorrectors)
            {
                if (isA<dynamicRefineFvMesh>(mesh))
                {
                    advector.surf().reconstruct();
                }

                mesh.update();

                if (mesh.changing())
                {
                    gh = (g & mesh.C()) - ghRef;

                    if (isA<dynamicRefineFvMesh>(mesh))
                    {
                        advector.surf().mapAlphaField();
                        alpha2 = 1.0 - alpha1;
                        alpha2.correctBoundaryConditions();
                        rho == alpha1*rho1 + alpha2*rho2;
                        rho.correctBoundaryConditions();
                        rho.oldTime() = rho;
                        alpha2.oldTime() = alpha2;
                    }

                    MRF.update();

                    if (correctPhi)
                    {
                        // Calculate absolute flux
                        // from the mapped surface velocity
                        phi = mesh.Sf() & Uf();

                        #include "correctPhi.H"

                        // Make the flux relative to the mesh motion
                        fvc::makeRelative(phi, U);

                        mixture.correct();
                    }

                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }
                }
            }

            #include "alphaControls.H"
            #include "alphaEqnSubCycle.H"

            mixture.correct();

            if (pimple.frozenFlow())
            {
                continue;
            }

            #include "gravityForce.H"
            #include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }

        runTime.write();
        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
