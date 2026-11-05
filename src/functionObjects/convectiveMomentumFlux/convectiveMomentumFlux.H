/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Class
    Foam::functionObjects::convectiveMomentumFlux

Description
    Face flux of momentum of the momentum convection term of the last time
    step, rhoPhi*U_f [kg m/s^2], with U_f the face velocity of the convection
    scheme of div(rhoPhi,U), weighted in time as in interIsoGoam:

        U_f = theta*interpolate(rhoPhi, U^{n+1})
            + (1 - theta)*interpolate(rhoPhi, U^n)

    theta is the PIMPLE entry convectionTheta, with the solver's defaults if
    absent: 1/2 for nOuterCorrectors > 1, else 0 with explicitConvection yes
    and 1 otherwise. A diagnostic: run it at write times (executeControl
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
            theta           <from PIMPLE>;
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

        //- Time weight (negative: from the PIMPLE dictionary)
        scalar theta_;

        //- Name of the result field
        word resultName_;


    // Private Member Functions

        //- Time weight of the convection term
        scalar theta() const;


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
