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

#include "interfaceUpwind.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::surfaceScalarField> Foam::interfaceUpwind<Type>::weights
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    const fvMesh& mesh = this->mesh();
    const volScalarField& alpha =
        mesh.template lookupObject<volScalarField>(alphaName_);

    tmp<surfaceScalarField> tw(tScheme_().weights(vf));
    surfaceScalarField& w = tw.ref();

    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();
    scalarField& wIn = w.primitiveFieldRef();
    const scalarField& phiIn = faceFlux_.primitiveField();

    const auto samePure = [this](const scalar aP, const scalar aN)
    {
        return (aP <= tol_ && aN <= tol_) || (aP >= 1 - tol_ && aN >= 1 - tol_);
    };

    forAll(wIn, facei)
    {
        if (!samePure(alpha[own[facei]], alpha[nei[facei]]))
        {
            wIn[facei] = pos0(phiIn[facei]);
        }
    }

    // Coupled patches
    auto& wbf = w.boundaryFieldRef();
    forAll(wbf, patchi)
    {
        const fvPatchScalarField& ap = alpha.boundaryField()[patchi];
        if (ap.coupled())
        {
            const scalarField aP(ap.patchInternalField());
            const scalarField aN(ap.patchNeighbourField());
            const scalarField& phip = faceFlux_.boundaryField()[patchi];
            scalarField& wp = wbf[patchi];
            forAll(wp, i)
            {
                if (!samePure(aP[i], aN[i]))
                {
                    wp[i] = pos0(phip[i]);
                }
            }
        }
    }

    return tw;
}


// ************************************************************************* //
