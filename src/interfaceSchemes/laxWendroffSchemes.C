/*---------------------------------------------------------------------------*\
    Copyright (C) 2026 Johan Roenby
    SPDX-License-Identifier: GPL-3.0-or-later

    Instantiation of the laxWendroff interpolation scheme (vector only)
\*---------------------------------------------------------------------------*/

#include "laxWendroff.H"
#include "fvMesh.H"

namespace Foam
{
    makeSurfaceInterpolationTypeScheme(laxWendroff, vector)
}

// ************************************************************************* //
