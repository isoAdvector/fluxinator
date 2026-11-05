/*---------------------------------------------------------------------------*\
    Instantiation of the laxWendroff interpolation scheme (vector only)
\*---------------------------------------------------------------------------*/

#include "laxWendroff.H"
#include "fvMesh.H"

namespace Foam
{
    makeSurfaceInterpolationTypeScheme(laxWendroff, vector)
}

// ************************************************************************* //
