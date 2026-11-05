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

Class
    Foam::basicGravityModel

Description
    Gravity force calculated as ghf*fvc::sngrad(rho).

SourceFiles
    basicGravityModel.C

\*---------------------------------------------------------------------------*/

#ifndef basicGravityModel_H
#define basicGravityModel_H

#include "typeInfo.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "autoPtr.H"
#include "uniformDimensionedFields.H"

#include "gravityModel.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class basicGravityModel Declaration
\*---------------------------------------------------------------------------*/

class basicGravityModel
:
public gravityModel
{

        //- gravity vector
        dimensionedVector g_;

        //- Href Field
        uniformDimensionedScalarField hRef_;

    // Private Member Functions

        //- Disallow copy construct
        basicGravityModel(const basicGravityModel&);

        //- Disallow default bitwise assignment
        void operator=(const basicGravityModel&);


public:

    //- Runtime type information
    TypeName("basic");


    // Constructors

        //- Construct from components
        basicGravityModel
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~basicGravityModel()
    {}

    virtual tmp<surfaceScalarField> force();


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
