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
    Foam::laxWendroff

Description
    Lax-Wendroff face value for an explicitly (start-of-step)
    evaluated convection term. The face value is advanced half a time step
    along the flow,
        vf_f = vf_f,s - lambda (dt/2) (u_f & grad(vf)_f),
    where vf_f,s is the face value of an uncorrected underlying scheme
    (linear: plain Lax-Wendroff; a limited scheme such as limitedLinearV:
    flux-limited Lax-Wendroff, Sweby 1984, with lambda the limiter of that
    scheme; interfaceUpwind: no correction on its upwind interface faces).
    lambda is the fraction of linear weights in the underlying weights,
    w = lambda w_lin + (1 - lambda) w_upwind. This cancels the
    anti-diffusion (dt/2)(u & grad)^2 of forward Euler time integration of
    the convection term (cf. Remmerswaal & Veldman 2022, J. Comput. Phys.,
    Lax-Wendroff flux interpolant).

    u_f: linear interpolate of U.oldTime() with the normal component replaced
    by phi/|Sf|; grad(vf)_f: linear interpolate of the Gauss gradient with the
    normal component replaced by the compact face-normal gradient.

    The correction is scaled by 1 - 2 w, with w the convectionNewTimeWeight
    that fluxinator and interIsoGoam register: the time-weighted convection
    has (1 - 2 w) times the anti-diffusion of forward Euler, so the
    correction vanishes for Crank-Nicolson convection (w 1/2). Without w
    (another solver) the scheme stops with an error.

    Usage, for explicit momentum convection in fluxinator, or in interIsoGoam
    with PIMPLE convectionNewTimeWeight 0 (convected field U.oldTime()):
    \verbatim
        div(rhoPhi,U)   Gauss laxWendroff interfaceUpwind limitedLinearV 1;
    \endverbatim
    The velocity field is U and the volumetric face flux phi, as in the
    solvers of OpenFOAM. A corrected underlying scheme is allowed: its
    correction is added to that of laxWendroff.
    fluxinator links the library; other solvers need libs (interfaceSchemes)
    in controlDict.

SourceFiles
    laxWendroff.C

\*---------------------------------------------------------------------------*/

#ifndef laxWendroff_H
#define laxWendroff_H

#include "surfaceInterpolationScheme.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
class laxWendroff
:
    public surfaceInterpolationScheme<Type>
{
    // Private Data

        //- Name of the velocity field (its oldTime() is used)
        const word UName_;

        //- Name of the volumetric face flux
        const word phiName_;

        //- Face flux of the underlying scheme (upwind direction)
        const surfaceScalarField& faceFlux_;

        //- Underlying scheme
        tmp<surfaceInterpolationScheme<Type>> tScheme_;


    // Private Member Functions

        //- The correction, given the weights w of the underlying scheme
        tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> correction
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf,
            const surfaceScalarField& w
        ) const;


public:

    //- Runtime type information
    TypeName("laxWendroff");


    // Constructors

        //- Construct from mesh and Istream
        laxWendroff(const fvMesh& mesh, Istream& is)
        :
            surfaceInterpolationScheme<Type>(mesh),
            UName_("U"),
            phiName_("phi"),
            faceFlux_(mesh.lookupObject<surfaceScalarField>(phiName_)),
            tScheme_
            (
                surfaceInterpolationScheme<Type>::New
                (
                    mesh,
                    mesh.lookupObject<surfaceScalarField>(phiName_),
                    is
                )
            )
        {}

        //- Construct from mesh, faceFlux and Istream
        laxWendroff
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& is
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            UName_("U"),
            phiName_("phi"),
            faceFlux_(faceFlux),
            tScheme_(surfaceInterpolationScheme<Type>::New(mesh, faceFlux, is))
        {}

        //- No copy construct
        laxWendroff(const laxWendroff&) = delete;

        //- No copy assignment
        void operator=(const laxWendroff&) = delete;


    // Member Functions

        //- Weights of the underlying scheme
        tmp<surfaceScalarField> weights
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const
        {
            return tScheme_().weights(vf);
        }

        //- Explicit Lax-Wendroff correction
        virtual bool corrected() const
        {
            return true;
        }

        //- -lambda (dt/2) (u_f & grad(vf)_f) on internal and coupled faces
        virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
        correction
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const;

        using surfaceInterpolationScheme<Type>::interpolate;

        //- Face values with the correction; the weights of the underlying
        //  scheme are computed once
        virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
        interpolate
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "laxWendroff.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
