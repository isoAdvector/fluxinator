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
    Foam::interfaceUpwind

Description
    Upwind interpolation on faces at the interface, any uncorrected
    interpolation scheme elsewhere. A face uses the underlying scheme only if
    both cells are pure cells of the same phase (alpha <= tol in both or
    alpha >= 1 - tol in both); on every other face (touching a mixed cell or
    separating the phases) the upwind value is used.

    With one velocity per cell, the velocity on the far side of the
    interface belongs to the other phase (the tangential velocity may jump
    at the interface), so the only face velocity consistent with the phase
    mass fluxes of a geometric VoF method is that of the cell the mass comes
    from, for the heavy and the light fluid alike. A central value would
    transfer momentum between the phases; a density-weighted value would be
    a downwind value for the light fluid.

    Usage, e.g. for the momentum convection of a VoF solver:
    \verbatim
        div(rhoPhi,U)   Gauss interfaceUpwind alpha.water 1e-6 limitedLinearV 1;
    \endverbatim
    fluxinator links the library; other solvers need libs (interfaceSchemes)
    in controlDict. Needs the face flux; the underlying scheme must not be a
    corrected scheme.

SourceFiles
    interfaceUpwind.C

\*---------------------------------------------------------------------------*/

#ifndef interfaceUpwind_H
#define interfaceUpwind_H

#include "surfaceInterpolationScheme.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
class interfaceUpwind
:
    public surfaceInterpolationScheme<Type>
{
    // Private Data

        //- Name of the volume fraction field
        const word alphaName_;

        //- Tolerance defining pure cells
        const scalar tol_;

        //- Face flux (for the upwind direction)
        const surfaceScalarField& faceFlux_;

        //- Underlying scheme
        tmp<surfaceInterpolationScheme<Type>> tScheme_;


    // Private Member Functions

        void checkUncorrected() const
        {
            if (tScheme_().corrected())
            {
                FatalErrorInFunction
                    << "interfaceUpwind needs an uncorrected underlying scheme"
                    << exit(FatalError);
            }
        }


public:

    //- Runtime type information
    TypeName("interfaceUpwind");


    // Constructors

        //- Construct from mesh and Istream: not available (needs the flux)
        interfaceUpwind(const fvMesh& mesh, Istream& is)
        :
            surfaceInterpolationScheme<Type>(mesh),
            alphaName_(is),
            tol_(readScalar(is)),
            faceFlux_
            (
                mesh.lookupObject<surfaceScalarField>(word(is))
            ),
            tScheme_(surfaceInterpolationScheme<Type>::New(mesh, faceFlux_, is))
        {
            checkUncorrected();
        }

        //- Construct from mesh, faceFlux and Istream
        interfaceUpwind
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& is
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            alphaName_(is),
            tol_(readScalar(is)),
            faceFlux_(faceFlux),
            tScheme_(surfaceInterpolationScheme<Type>::New(mesh, faceFlux, is))
        {
            checkUncorrected();
        }

        //- No copy construct
        interfaceUpwind(const interfaceUpwind&) = delete;

        //- No copy assignment
        void operator=(const interfaceUpwind&) = delete;


    // Member Functions

        //- Owner weights: upwind at interface faces, underlying scheme elsewhere
        tmp<surfaceScalarField> weights
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "interfaceUpwind.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
