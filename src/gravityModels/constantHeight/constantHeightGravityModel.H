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
    Foam::constantHeightGravityModel

Description
    Gravity force with the height of a fixed horizontal interface on every
    face, gh_f = g & interfacePoint - ghRef: the exact gravity term for a flat
    interface at rest, a reference for testing the other models.

Usage
    In system/fvSchemes:
    \verbatim
    interfaceForces
    {
        gravity
        {
            model           constantHeight;
            interfacePoint  (0 0 0.5);  // a point on the interface
        }
    }
    \endverbatim

SourceFiles
    constantHeightGravityModel.C

\*---------------------------------------------------------------------------*/

#ifndef constantHeightGravityModel_H
#define constantHeightGravityModel_H

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
                Class constantHeightGravityModel Declaration
\*---------------------------------------------------------------------------*/

class constantHeightGravityModel
:
    public gravityModel
{
    // Private Data

        //- A point on the interface
        const point interfacePoint_;

        //- g*h of the interface, used on every face
        const dimensionedScalar ghInterface_;


public:

    //- Runtime type information
    TypeName("constantHeight");


    // Constructors

        //- Construct from components
        constantHeightGravityModel
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~constantHeightGravityModel() = default;


    // Member Functions

        //- The face gravity force
        virtual tmp<surfaceScalarField> force();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
