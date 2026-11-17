/*---------------------------------------------------------------------------*\
            Copyright (c) 2017-2019, German Aerospace Center (DLR)
    Modified code Copyright (C) 2026 Johan Roenby
-------------------------------------------------------------------------------
License
    This file is part of the VoFLibrary source code library, which is an
    unofficial extension to OpenFOAM.

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
    Foam::markInterfaceRegion

Description
    Marks the interface cells and a given number of layers of point neighbours
    around them, with the layer number of each cell.

Author
    Henning Scheufler, DLR, all rights reserved.

SourceFiles
    markInterfaceRegion.C

\*---------------------------------------------------------------------------*/

#ifndef markInterfaceRegion_H
#define markInterfaceRegion_H

#include "fvMesh.H"
#include "indirectPrimitivePatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class markInterfaceRegion Declaration
\*---------------------------------------------------------------------------*/

class markInterfaceRegion
{
    // Private Data

        //- The mesh
        const fvMesh& mesh_;

        //- Points on coupled boundaries, which are synchronised
        labelList coupledBoundaryPoints_;


    // Private Member Functions

        //- Patch of all coupled faces
        autoPtr<indirectPrimitivePatch> coupledFacesPatch() const;


public:

    // Constructors

        //- Construct from mesh
        explicit markInterfaceRegion(const fvMesh& mesh);


    // Member Functions

        //- Mark the interface cells (layer 0) and nLayers layers of point
        //  neighbours around them. Layer k consists of the cells that share a
        //  point with layer k - 1 and are not in a lower layer. On return,
        //  nearInterface is true in all of them, and cellLayer holds the layer
        //  of each cell, -1 outside.
        void markCellsNearSurf
        (
            const boolList& interfaceCells,
            const label nLayers,
            boolList& nearInterface,
            labelList& cellLayer
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
