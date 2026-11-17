/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Johan Roenby
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
    Foam::functionObjects::convectiveMomentumFlux

Description
    Face flux of momentum of the momentum convection term of the last time
    step, rhoPhi*U_f [kg m/s^2], with U_f the face velocity of the convection
    scheme of div(rhoPhi,U), weighted in time as in fluxinator and interIsoGoam:

        U_f = w*interpolate(rhoPhi, U^{n+1})
            + (1 - w)*interpolate(rhoPhi, U^n)

    w is the convectionNewTimeWeight that fluxinator and interIsoGoam
    register. With other solvers it is the entry convectionNewTimeWeight of
    this function object, default 1 (implicit convection, as in the solvers
    of OpenFOAM). A diagnostic: run it at write times (executeControl
    writeTime), when U^n is still the old-time velocity of the step.

    Usage, in system/controlDict functions:
    \verbatim
        momentumFlux
        {
            type            convectiveMomentumFlux;
            libs            (interfaceFunctionObjects);
            executeControl  writeTime;
            writeControl    writeTime;

            // Optional entries (defaults)
            U               U;
            rhoPhi          rhoPhi;
            scheme          div(rhoPhi,U);
            convectionNewTimeWeight 1;  // only for other solvers
            result          convectiveMomentumFlux;
        }
    \endverbatim

SourceFiles
    convectiveMomentumFlux.C

\*---------------------------------------------------------------------------*/

#ifndef convectiveMomentumFlux_H
#define convectiveMomentumFlux_H

#include "fvMeshFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

class convectiveMomentumFlux
:
    public fvMeshFunctionObject
{
    // Private Data

        //- Name of the velocity field
        word UName_;

        //- Name of the mass flux field
        word rhoPhiName_;

        //- Name of the convection scheme in fvSchemes divSchemes
        word schemeName_;

        //- Entry convectionNewTimeWeight, for solvers that do not register it
        scalar weight_;

        //- Whether the entry was given
        bool weightGiven_;

        //- Name of the result field
        word resultName_;


    // Private Member Functions

        //- Time weight of the convection term
        scalar weight() const;


public:

    //- Runtime type information
    TypeName("convectiveMomentumFlux");


    // Constructors

        //- Construct from name, Time and dictionary
        convectiveMomentumFlux
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        convectiveMomentumFlux(const convectiveMomentumFlux&) = delete;

        //- No copy assignment
        void operator=(const convectiveMomentumFlux&) = delete;


    //- Destructor
    virtual ~convectiveMomentumFlux() = default;


    // Member Functions

        //- Read the settings
        virtual bool read(const dictionary& dict);

        //- Calculate the flux
        virtual bool execute();

        //- Write the flux
        virtual bool write();
};


} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
