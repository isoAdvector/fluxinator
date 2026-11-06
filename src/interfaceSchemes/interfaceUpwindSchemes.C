/*---------------------------------------------------------------------------*\
    Copyright (C) 2026 Johan Roenby
    SPDX-License-Identifier: GPL-3.0-or-later

    Instantiation of the interfaceUpwind interpolation scheme for all types
\*---------------------------------------------------------------------------*/

#include "interfaceUpwind.H"
#include "fvMesh.H"

namespace Foam
{
    makeSurfaceInterpolationScheme(interfaceUpwind)
}

// ************************************************************************* //
