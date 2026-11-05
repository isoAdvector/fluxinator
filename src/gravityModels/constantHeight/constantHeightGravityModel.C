/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2025 Johan Roenby
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

\*---------------------------------------------------------------------------*/

#include "constantHeightGravityModel.H"
#include "gravityMeshObject.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(constantHeightGravityModel, 0);
    addToRunTimeSelectionTable
    (
        gravityModel,
        constantHeightGravityModel,
        components
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::constantHeightGravityModel::constantHeightGravityModel
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    gravityModel
    (
        typeName,
        dict,
        mesh
    ),
    ghf_
    (
        IOobject
        (
            "ghf",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("ghf_", dimAcceleration*dimLength, 0.0)
    )
{

    // Reading surface point
    IOdictionary setAlphaFieldDict
    (
        IOobject
        (
            "setAlphaFieldDict",
            mesh_.time().system(),
            mesh_.time(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    );
    point interfacePoint = (setAlphaFieldDict.get<vector>("origin"));

    // Reading gravity
    const uniformDimensionedVectorField& g =
        meshObjects::gravity::New(mesh_.time());

    // Reading reference height or setting to zero
    uniformDimensionedScalarField hRef
    (
        IOobject
        (
            "hRef",
            mesh.time().constant(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        dimensionedScalar(dimLength, Zero)
    );

    dimensionedScalar ghRef
    (
        mag(g.value()) > SMALL
        ? g & (cmptMag(g.value())/mag(g.value()))*hRef
        : dimensionedScalar("ghRef", g.dimensions()*dimLength, 0)
    );

    // Setting internal field
    ghf_.primitiveFieldRef() = (g.value() & interfacePoint) - ghRef.value();

    // Setting boundary field
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    forAll(patches, patchi)
    {
        ghf_.boundaryFieldRef()[patchi] =
            (g.value() & interfacePoint) - ghRef.value();
    }
}


// * * * * * * * * * * * * * * Public Access Member Functions  * * * * * * * //

Foam::tmp<Foam::surfaceScalarField> Foam::constantHeightGravityModel::force()
{
    const volScalarField& rho = mesh_.lookupObject<volScalarField>("rho");
    return -ghf_*fvc::snGrad(rho);
}

// ************************************************************************* //
