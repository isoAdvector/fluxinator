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

Description
    The reference geopotential ghRef = g & (|g|/|g|) hRef of the VoF solvers
    of OpenFOAM (gh.H): the reference height hRef (constant/hRef) measured
    along the components of gravity, zero without gravity. gh = g & x - ghRef
    is then the geopotential relative to it. Used by fluxinator, its gravity
    models and setHydrostaticPressure, so that they all use the same value.

\*---------------------------------------------------------------------------*/

#ifndef ghRef_H
#define ghRef_H

#include "dimensionedScalar.H"
#include "dimensionedVector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- The reference geopotential for gravity g and reference height hRef
inline scalar ghRef(const vector& g, const scalar hRef)
{
    return mag(g) > SMALL ? g & (cmptMag(g)/mag(g))*hRef : 0;
}

//- The reference geopotential for gravity g and reference height hRef
inline dimensionedScalar ghRef
(
    const dimensionedVector& g,
    const dimensionedScalar& hRef
)
{
    return dimensionedScalar
    (
        "ghRef",
        g.dimensions()*hRef.dimensions(),
        ghRef(g.value(), hRef.value())
    );
}

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
