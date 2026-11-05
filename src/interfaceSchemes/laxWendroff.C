/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#include "laxWendroff.H"
#include "surfaceFields.H"
#include "linear.H"
#include "gaussGrad.H"
#include "snGradScheme.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::laxWendroff<Type>::correction
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    typedef typename outerProduct<vector, Type>::type GradType;

    const fvMesh& mesh = this->mesh();
    const volVectorField& U =
        mesh.template lookupObject<volVectorField>(UName_).oldTime();
    const surfaceScalarField& phi =
        mesh.template lookupObject<surfaceScalarField>(phiName_);
    const surfaceVectorField nf(mesh.Sf()/mesh.magSf());

    // Face velocity: tangential part interpolated, normal part from the flux
    surfaceVectorField uf(linearInterpolate(U));
    uf += nf*(phi/mesh.magSf() - (nf & uf));

    // Face gradient: tangential part interpolated, normal part compact
    GeometricField<GradType, fvsPatchField, surfaceMesh> gradf
    (
        linearInterpolate(fv::gaussGrad<Type>::gradf(linearInterpolate(vf), "gradLW"))
    );
    gradf +=
        nf
       *(
            fv::snGradScheme<Type>::snGrad(vf, mesh.nonOrthDeltaCoeffs())
          - (nf & gradf)
        );

    // Time-weighted convection theta*div(U^{n+1}) + (1 - theta)*div(U^n):
    // the anti-diffusion is (1 - 2 theta) times that of forward Euler, and
    // the solver registers laxWendroffFactor = 1 - 2 theta (1 if absent)
    const scalar lwFactor =
        mesh.template foundObject<uniformDimensionedScalarField>
        (
            "laxWendroffFactor"
        )
      ? mesh.template lookupObject<uniformDimensionedScalarField>
        (
            "laxWendroffFactor"
        ).value()
      : 1;

    const dimensionedScalar dt(lwFactor*mesh.time().deltaT());

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tcorr
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                "laxWendroff::correction(" + vf.name() + ')',
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            -0.5*dt*(uf & gradf)
        )
    );
    GeometricField<Type, fvsPatchField, surfaceMesh>& corr = tcorr.ref();

    // Flux-limited variant: scale by the fraction lambda of the linear
    // weights in the underlying weights, w = lambda w_lin + (1 - lambda)
    // w_upwind (the limiter of a limited scheme, 0 on the upwind faces of
    // interfaceUpwind, 1 for linear)
    {
        const surfaceScalarField w(tScheme_().weights(vf));
        const surfaceScalarField& wLin = mesh.weights();
        const auto lambda = [](scalar w, scalar wl, scalar fl)
        {
            const scalar wu = pos0(fl);
            return mag(wl - wu) > SMALL ? (w - wu)/(wl - wu) : 1.0;
        };
        Field<Type>& cIn = corr.primitiveFieldRef();
        forAll(cIn, facei)
        {
            cIn[facei] *=
                lambda(w[facei], wLin[facei], faceFlux_[facei]);
        }
        auto& cbfL = corr.boundaryFieldRef();
        forAll(cbfL, patchi)
        {
            if (vf.boundaryField()[patchi].coupled())
            {
                const scalarField& wp = w.boundaryField()[patchi];
                const scalarField& wlp = wLin.boundaryField()[patchi];
                const scalarField& fp = faceFlux_.boundaryField()[patchi];
                forAll(cbfL[patchi], i)
                {
                    cbfL[patchi][i] *= lambda(wp[i], wlp[i], fp[i]);
                }
            }
        }
    }

    // No correction on non-coupled boundaries (boundary values are used)
    auto& cbf = corr.boundaryFieldRef();
    forAll(cbf, patchi)
    {
        if (!vf.boundaryField()[patchi].coupled())
        {
            cbf[patchi] = Zero;
        }
    }

    // Correction of a corrected underlying scheme
    if (tScheme_().corrected())
    {
        corr += tScheme_().correction(vf);
    }

    return tcorr;
}


// ************************************************************************* //
