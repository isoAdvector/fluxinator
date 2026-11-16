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
#include "fvSchemes.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::gravityModel>
Foam::gravityModel::New
(
    const fvMesh& mesh,
    const word& defaultModel
)
{
    // interfaceForces.gravity in system/fvSchemes, read as OpenFOAM reads
    // wallDist there; without it the solver's default model
    const dictionary dict
    (
        static_cast<const fvSchemes&>(mesh)
            .subOrEmptyDict("interfaceForces")
            .subOrEmptyDict("gravity")
    );

    const word modelType(dict.getOrDefault<word>("model", defaultModel));

    Info<< "Selecting gravity model " << modelType
        << (dict.found("model") ? " (" : " (default, no ")
        << "interfaceForces.gravity.model in system/fvSchemes)" << endl;

    auto* ctorPtr = componentsConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "model",
            modelType,
            *componentsConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<gravityModel>(ctorPtr(dict, mesh));
}


// ************************************************************************* //
