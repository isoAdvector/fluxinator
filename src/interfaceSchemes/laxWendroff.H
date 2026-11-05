/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

    If the mesh holds a uniformDimensionedScalarField laxWendroffFactor
    (interIsoGoam: 1 - 2 convectionTheta), the correction is scaled by
    it: the theta scheme has (1 - 2 theta) times the anti-diffusion of
    forward Euler, so the correction vanishes for the time-centred scheme.

    Usage, for explicit momentum convection in interIsoGoam
    (PIMPLE explicitConvection yes, convected field U.oldTime()):
    \verbatim
        div(rhoPhi,U)   Gauss laxWendroff U phi limitedLinearV 1;
    \endverbatim
    with libs (interfaceSchemes) in controlDict.

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

        //- A corrected underlying scheme is allowed: its correction is
        //  added to that of laxWendroff (the limiter fraction
        //  is taken from its weights)
        void checkUncorrected() const
        {}

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
            UName_(is),
            phiName_(is),
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
        {
            checkUncorrected();
        }

        //- Construct from mesh, faceFlux and Istream
        laxWendroff
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& is
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            UName_(is),
            phiName_(is),
            faceFlux_(faceFlux),
            tScheme_(surfaceInterpolationScheme<Type>::New(mesh, faceFlux, is))
        {
            checkUncorrected();
        }

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
