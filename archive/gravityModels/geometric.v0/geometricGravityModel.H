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
    Foam::geometricGravityModel

Description
    Gravity force calculated based on reconstructed interface centres
    interpolated to mesh faces.

SourceFiles
    geometricGravityModel.C

\*---------------------------------------------------------------------------*/

#ifndef geometricGravityModel_H
#define geometricGravityModel_H

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
                Class geometricGravityModel Declaration
\*---------------------------------------------------------------------------*/

class geometricGravityModel
:
    public gravityModel
{

        // g*h at faces based on geometrically reconstructed interface height
        surfaceScalarField ghf_;

        //- Height reference field
        uniformDimensionedScalarField hRef_;

    // Private Member Functions

        //- Disallow copy construct
        geometricGravityModel(const geometricGravityModel&);

        //- Disallow default bitwise assignment
        void operator=(const geometricGravityModel&);


        protected:

        //- Calculate the interface height
        void calculateGhf();


public:

    //- Runtime type information
    TypeName("geometric");


    // Constructors

        //- Construct from components
        geometricGravityModel
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~geometricGravityModel()
    {}

    virtual tmp<surfaceScalarField> force();


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
