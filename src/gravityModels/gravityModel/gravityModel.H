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

Class
    Foam::gravityModel

Description
    Abstract base class for models of the gravity force in interfacial
    solvers like interFoam and interIsoFoam: the face force -gh_f snGrad(rho)
    that remains when the hydrostatic part is absorbed in p_rgh, with
    gh = g & x - ghRef.

SourceFiles
    gravityModel.C
    gravityModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef gravityModel_H
#define gravityModel_H

#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "autoPtr.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "uniformDimensionedFields.H"
#include "ghRef.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class gravityModel Declaration
\*---------------------------------------------------------------------------*/

class gravityModel
{
protected:

    // Protected Data

        const fvMesh& mesh_;

        //- Reference height of the solvers (constant/hRef, default 0)
        const uniformDimensionedScalarField hRef_;


    // Protected Member Functions

        //- The reference geopotential of the solvers (ghRef.H)
        dimensionedScalar ghRef() const;


public:

    //- Runtime type information
    TypeName("gravityModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            gravityModel,
            components,
            (
                const dictionary& dict,
                const fvMesh& mesh
            ),
            (dict, mesh)
        );


    // Selectors

        //- Return a reference to the selected gravity model: the model of
        //  interfaceForces.gravity in system/fvSchemes, or defaultModel
        //  without that entry
        static autoPtr<gravityModel> New
        (
            const fvMesh& mesh,
            const word& defaultModel = "basic"
        );


    // Constructors

        //- Construct from components
        gravityModel
        (
            const word& type,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- No copy construct
        gravityModel(const gravityModel&) = delete;

        //- No copy assignment
        void operator=(const gravityModel&) = delete;


    //- Destructor
    virtual ~gravityModel() = default;


    // Member Functions

        //- The face gravity force
        virtual tmp<surfaceScalarField> force() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
