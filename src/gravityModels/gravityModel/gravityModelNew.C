/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2025-2026 Johan Roenby
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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::gravityModel>
Foam::gravityModel::New
(
    const fvMesh& mesh,
    const word& defaultModel
)
{
    // Create an IOdictionary with default values
    IOdictionary gravityProperties
    (
        IOobject
        (
            "gravityProperties",
            mesh.time().constant(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false  // Not registered
        )
    );

    // Add default value if the file wasn't found
    if (!gravityProperties.headerOk())
    {
        gravityProperties.add("gravityModel", defaultModel);
        Info<< "No gravityProperties file found. Using default values." << endl;
    }
    else
    {
        Info<< "Reading gravityProperties from constant/gravityProperties"
            << endl;
    }

    word gravityModelTypeName = gravityProperties.get<word>("gravityModel");
    Info<< "Selecting gravityModel: " << gravityModelTypeName << endl;

    auto* ctorPtr = componentsConstructorTable(gravityModelTypeName);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            gravityProperties,
            "gravityModel",
            gravityModelTypeName,
            *componentsConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<gravityModel>(ctorPtr(gravityProperties, mesh));
}


// ************************************************************************* //
