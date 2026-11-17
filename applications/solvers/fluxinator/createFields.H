// Modified code Copyright (C) 2026 Johan Roenby
// SPDX-License-Identifier: GPL-3.0-or-later

#include "createRDeltaT.H"

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"


Info<< "Reading transportProperties\n" << endl;
immiscibleIncompressibleTwoPhaseMixture mixture(U, phi);

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());

const dimensionedScalar& rho1 = mixture.rho1();
const dimensionedScalar& rho2 = mixture.rho2();


// Need to store rho for ddt(rho, U)
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    alpha1*rho1 + alpha2*rho2
);
rho.oldTime();


// Mass flux
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
);


typedef incompressibleInterPhaseTransportModel
    <
        immiscibleIncompressibleTwoPhaseMixture
    > transportModelType;

autoPtr<transportModelType> turbulence
(
    new transportModelType
    (
        rho, U, phi, rhoPhi, mixture
    )
);


#include "readGravitationalAcceleration.H"
#include "readhRef.H"

// The geopotential at the cell centres, for p = p_rgh + rho*gh, as gh.H of
// OpenFOAM; the gravity models compute the one at the faces
Info<< "Calculating field g.h\n" << endl;
const dimensionedScalar ghRef(Foam::ghRef(g, hRef));
const auto oldConsistency = FieldBase::localBoundaryConsistency(0);
volScalarField gh("gh", (g & mesh.C()) - ghRef);
FieldBase::localBoundaryConsistency(oldConsistency);


volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rho*gh
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);

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

mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha1.name());

#include "createMRF.H"
#include "createFvOptions.H"
#include "createPorosity.H"

isoAdvector advector(alpha1, phi, U);

// Reconstructions of the initial interface before the first time step. An
// iterative scheme such as plicRDF starts each reconstruction from the normals
// of the previous one, so they refine the interface that the first advection
// step starts from (solver dictionary of alpha, default in solverDefaults.H).
{
    const label nInitialReconstructions =
        mesh.solverDict(alpha1.name()).getOrDefault<label>
        (
            "nInitialReconstructions",
            defaultNInitialReconstructions
        );

    Info<< "Initial interface: nInitialReconstructions "
        << nInitialReconstructions << endl;

    for (label i = 0; i < nInitialReconstructions; ++i)
    {
        advector.surf().reconstruct(true);
    }
}

// Gravity model: interfaceForces.gravity in system/fvSchemes, without it the
// default of solverDefaults.H
autoPtr<gravityModel> gravity = gravityModel::New(mesh, defaultGravityModel);

// Run-time controls, read once per time step (readControls.H). The convection
// time weight is also registered, for the laxWendroff scheme and the
// convectiveMomentumFlux function object.
scalar convectionWeight = -1;
scalar gravityWeight = -1;
scalar outerWeight = -1;
bool ddtCorr = defaultDdtCorr;
bool convectionSchemeChecked = false;
uniformDimensionedScalarField convectionNewTimeWeight
(
    IOobject
    (
        "convectionNewTimeWeight",
        runTime.constant(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        IOobject::REGISTER
    ),
    dimensionedScalar(dimless, Zero)
);
#include "readControls.H"

// The gravity force of the latest outer corrector and that of the last outer
// corrector of the previous time step, for the time weighting (gravityForce.H).
// Only needed with gravityNewTimeWeight != 1 (outer correctors).
autoPtr<surfaceScalarField> gravityFLatestPtr;
autoPtr<surfaceScalarField> gravityFOldPtr;
if (gravityWeight != 1)
{
    gravityFLatestPtr.reset
    (
        new surfaceScalarField("gravityFLatest", gravity->force())
    );
    gravityFOldPtr.reset
    (
        new surfaceScalarField("gravityFOld", gravityFLatestPtr())
    );
}

