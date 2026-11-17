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
        div(rhoPhi,U)   Gauss interfaceUpwind limitedLinearV 1;
    \endverbatim
    The volume fraction field is alpha1 (alpha.<first phase> of
    transportProperties) and tol is surfCellTol of its solver dictionary
    (default 1e-8), so that the interface cells are those of isoAdvector.
    The phase makes no difference: the test is the same for alpha and
    1 - alpha.
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
#include "IOdictionary.H"

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

        //- alpha.<first phase> of constant/transportProperties, the alpha1
        //  of a two-phase mixture
        static word alpha1Name(const fvMesh& mesh)
        {
            const word dictName("transportProperties");
            const IOdictionary* dictPtr =
                mesh.findObject<IOdictionary>(dictName);
            if (dictPtr)
            {
                return IOobject::groupName
                (
                    "alpha",
                    dictPtr->get<wordList>("phases")[0]
                );
            }

            const IOdictionary dict
            (
                IOobject
                (
                    dictName,
                    mesh.time().constant(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                )
            );
            return IOobject::groupName("alpha", dict.get<wordList>("phases")[0]);
        }

        //- surfCellTol of the solver dictionary of alphaName, which defines
        //  the interface cells of isoAdvector (default 1e-8, as there)
        static scalar surfCellTol(const fvMesh& mesh, const word& alphaName)
        {
            const scalar defaultTol = 1e-8;
            const dictionary* dictPtr = mesh.solversDict().findDict(alphaName);
            return
                dictPtr
              ? dictPtr->getOrDefault<scalar>("surfCellTol", defaultTol)
              : defaultTol;
        }

        //- Fatal error for a corrected underlying scheme
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

        //- Construct from mesh and Istream, which names the face flux first
        interfaceUpwind(const fvMesh& mesh, Istream& is)
        :
            surfaceInterpolationScheme<Type>(mesh),
            alphaName_(alpha1Name(mesh)),
            tol_(surfCellTol(mesh, alphaName_)),
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
            alphaName_(alpha1Name(mesh)),
            tol_(surfCellTol(mesh, alphaName_)),
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
