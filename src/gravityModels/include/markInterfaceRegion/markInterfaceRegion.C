/*---------------------------------------------------------------------------*\
            Copyright (c) 2017-2019, German Aerospace Center (DLR)
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

\*---------------------------------------------------------------------------*/

//#include "dummyTransform.H"
#include "emptyPolyPatch.H"
#include "markInterfaceRegion.H"
#include "coupledPolyPatch.H"
#include "syncTools.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::autoPtr<Foam::indirectPrimitivePatch>
Foam::markInterfaceRegion::coupledFacesPatch() const
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    label nCoupled = 0;

    for (const polyPatch& pp : patches)
    {
        if (isA<coupledPolyPatch>(pp))
        {
            nCoupled += pp.size();
        }
    }
    labelList nCoupledFaces(nCoupled);
    nCoupled = 0;

    for (const polyPatch& pp : patches)
    {
        if (isA<coupledPolyPatch>(pp))
        {
            label facei = pp.start();

            forAll(pp, i)
            {
                nCoupledFaces[nCoupled++] = facei++;
            }
        }
    }

    return autoPtr<indirectPrimitivePatch>::New
    (
        IndirectList<face>
        (
            mesh_.faces(),
            nCoupledFaces
        ),
        mesh_.points()
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::markInterfaceRegion::markInterfaceRegion(const fvMesh& mesh)
:
    mesh_(mesh),
    coupledBoundaryPoints_(coupledFacesPatch()().meshPoints())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::markInterfaceRegion::markCellsNearSurf
(
    const boolList& interfaceCells,
    const label nLayers,
    boolList& nearInterface,
    labelList& cellLayer
)
{
    const labelListList& cellPoints = mesh_.cellPoints();
    const labelListList& pointCells = mesh_.pointCells();

    nearInterface.resize(mesh_.nCells());
    nearInterface = false;
    cellLayer.resize(mesh_.nCells());

    boolList alreadyMarkedPoint(mesh_.nPoints(), false);

    // Points on coupled boundaries next to the previous layer, synchronised
    // across processors
    Map<bool> syncMap;

    for (label layer = 0; layer <= nLayers; ++layer)
    {
        if (layer > 0)
        {
            // Mark the cells at the synchronised points first
            for (const label pointi : coupledBoundaryPoints_)
            {
                for (const label celli : pointCells[pointi])
                {
                    if (cellLayer[celli] == layer - 1)
                    {
                        syncMap.insert(pointi, true);
                        break;
                    }
                }
            }

            syncTools::syncPointMap(mesh_, syncMap, orEqOp<bool>());

            forAllConstIters(syncMap, iter)
            {
                const label pointi = iter.key();

                if (!alreadyMarkedPoint[pointi])
                {
                    for (const label celli : pointCells[pointi])
                    {
                        if (cellLayer[celli] == -1)
                        {
                            cellLayer[celli] = layer;
                            nearInterface[celli] = true;
                        }
                    }
                }

                alreadyMarkedPoint[pointi] = true;
            }
        }

        forAll(cellLayer, celli)
        {
            if (layer == 0)
            {
                if (interfaceCells[celli])
                {
                    cellLayer[celli] = 0;
                    nearInterface[celli] = true;
                }
                else
                {
                    cellLayer[celli] = -1;
                }
            }
            else if (cellLayer[celli] == layer - 1)
            {
                for (const label pointi : cellPoints[celli])
                {
                    if (!alreadyMarkedPoint[pointi])
                    {
                        for (const label pCelli : pointCells[pointi])
                        {
                            if (cellLayer[pCelli] == -1)
                            {
                                cellLayer[pCelli] = layer;
                                nearInterface[pCelli] = true;
                            }
                        }
                    }
                    alreadyMarkedPoint[pointi] = true;
                }
            }
        }
    }
}


// ************************************************************************* //
