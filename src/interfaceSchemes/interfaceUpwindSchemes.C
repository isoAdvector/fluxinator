/*---------------------------------------------------------------------------*\
    Instantiation of the interfaceUpwind interpolation scheme for all types
\*---------------------------------------------------------------------------*/

#include "interfaceUpwind.H"
#include "fvMesh.H"

namespace Foam
{
    makeSurfaceInterpolationScheme(interfaceUpwind)
}

// ************************************************************************* //
