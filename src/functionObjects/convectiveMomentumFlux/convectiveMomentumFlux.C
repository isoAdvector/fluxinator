/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#include "convectiveMomentumFlux.H"
#include "convectionScheme.H"
#include "uniformDimensionedFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(convectiveMomentumFlux, 0);
    addToRunTimeSelectionTable
    (
        functionObject,
        convectiveMomentumFlux,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::functionObjects::convectiveMomentumFlux::theta() const
{
    if (theta_ >= 0)
    {
        return theta_;
    }

    // The defaults of interIsoGoam (timeWeights.H)
    const dictionary& pimple = mesh_.solutionDict().subDict("PIMPLE");
    const bool centredOuter =
        pimple.getOrDefault<label>("nOuterCorrectors", 1) > 1;
    const bool explicitConvection =
        pimple.getOrDefault<bool>("explicitConvection", false);
    const scalar theta =
        pimple.getOrDefault<scalar>
        (
            "convectionTheta",
            centredOuter ? 0.5 : (explicitConvection ? 0 : 1)
        );

    // The solver registers 1 - 2 theta for the laxWendroff scheme
    const auto* lwPtr =
        mesh_.findObject<uniformDimensionedScalarField>("laxWendroffFactor");
    if (lwPtr && mag(lwPtr->value() - (1 - 2*theta)) > SMALL)
    {
        FatalErrorInFunction
            << "convectionTheta " << theta << " from the PIMPLE dictionary"
            << " differs from the solver's (laxWendroffFactor "
            << lwPtr->value() << "); set theta" << exit(FatalError);
    }

    return theta;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::convectiveMomentumFlux::convectiveMomentumFlux
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::convectiveMomentumFlux::read
(
    const dictionary& dict
)
{
    fvMeshFunctionObject::read(dict);

    UName_ = dict.getOrDefault<word>("U", "U");
    rhoPhiName_ = dict.getOrDefault<word>("rhoPhi", "rhoPhi");
    schemeName_ = dict.getOrDefault<word>("scheme", "div(rhoPhi,U)");
    theta_ = dict.getOrDefault<scalar>("theta", -1);
    resultName_ = dict.getOrDefault<word>("result", typeName);

    return true;
}


bool Foam::functionObjects::convectiveMomentumFlux::execute()
{
    const volVectorField& U = lookupObject<volVectorField>(UName_);
    const surfaceScalarField& rhoPhi =
        lookupObject<surfaceScalarField>(rhoPhiName_);

    tmp<fv::convectionScheme<vector>> tcs
    (
        fv::convectionScheme<vector>::New
        (
            mesh_, rhoPhi, mesh_.divScheme(schemeName_)
        )
    );

    const scalar t = theta();

    Log << type() << " " << name() << ": theta " << t << endl;

    return store
    (
        resultName_,
        rhoPhi
       *(
            t*tcs().interpolate(rhoPhi, U)
          + (1 - t)*tcs().interpolate(rhoPhi, U.oldTime())
        )
    );
}


bool Foam::functionObjects::convectiveMomentumFlux::write()
{
    return writeObject(resultName_);
}


// ************************************************************************* //
