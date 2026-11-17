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
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const surfaceScalarField& w
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
        linearInterpolate
        (
            fv::gaussGrad<Type>::gradf(linearInterpolate(vf), "gradLW")
        )
    );
    gradf +=
        nf
       *(
            fv::snGradScheme<Type>::snGrad(vf, mesh.nonOrthDeltaCoeffs())
          - (nf & gradf)
        );

    // Time-weighted convection w*div(U^{n+1}) + (1 - w)*div(U^n): the
    // anti-diffusion is (1 - 2 w) times that of forward Euler. The solver
    // registers w as convectionNewTimeWeight.
    const auto* weightPtr =
        mesh.template findObject<uniformDimensionedScalarField>
        (
            "convectionNewTimeWeight"
        );
    if (!weightPtr)
    {
        FatalErrorInFunction
            << "laxWendroff needs the time weight of the convection term,"
            << " convectionNewTimeWeight, which fluxinator and interIsoGoam"
            << " register" << exit(FatalError);
    }
    const scalar lwFactor = 1 - 2*weightPtr->value();

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


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::laxWendroff<Type>::correction
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    return correction(vf, tScheme_().weights(vf));
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::laxWendroff<Type>::interpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    // As surfaceInterpolationScheme::interpolate, with the weights of the
    // underlying scheme computed once for the face values and the correction
    const tmp<surfaceScalarField> tw(tScheme_().weights(vf));

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tsf
    (
        surfaceInterpolationScheme<Type>::interpolate
        (
            vf,
            tmp<surfaceScalarField>(tw())
        )
    );
    tsf.ref() += correction(vf, tw());

    return tsf;
}


// ************************************************************************* //
