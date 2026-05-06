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

#include "gravityModel.H"
#include "reconstructionSchemes.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(gravityModel, 0);
    defineRunTimeSelectionTable(gravityModel, components);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::gravityModel::gravityModel
(
    const word& type,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    mesh_(mesh)/*,
    accf_
    (
        IOobject
        (
            "ghf",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("accf_", dimAcceleration*dimLength, 0.0)
    ),
    acc_
    (
        IOobject
        (
            "gh",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("acc_", dimAcceleration*dimLength, 0.0)
    )
*/
{
}


// * * * * * * * * * * * * * * Public Access Member Functions  * * * * * * * //

/*
void Foam::gravityModel::calculateAcc()
{
    notImplemented("bool Foam::gravityModel::calculateAcc()");;
}
*/


// ************************************************************************* //
