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
    gravityModel(typeName, dict, mesh),
    interfacePoint_(dict.get<point>("interfacePoint")),
    ghInterface_
    (
        "ghInterface",
        dimAcceleration*dimLength,
        (meshObjects::gravity::New(mesh.time()).value() & interfacePoint_)
      - ghRef().value()
    )
{
    Info<< "constantHeight gravity model: interfacePoint " << interfacePoint_
        << endl;
}


// * * * * * * * * * * * * * * Public Access Member Functions  * * * * * * * //

Foam::tmp<Foam::surfaceScalarField> Foam::constantHeightGravityModel::force()
{
    const volScalarField& rho = mesh_.lookupObject<volScalarField>("rho");
    return -ghInterface_*fvc::snGrad(rho);
}

// ************************************************************************* //
