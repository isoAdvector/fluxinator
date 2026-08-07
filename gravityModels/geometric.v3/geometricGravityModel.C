/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2025 Johan Roenby
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

#include "geometricGravityModel.H"
#include "gravityMeshObject.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"
#include "reconstructionSchemes.H"
#include "cutFaceAdvect.H"
#include "processorPolyPatch.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(geometricGravityModel, 0);
    addToRunTimeSelectionTable(gravityModel, geometricGravityModel, components);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::geometricGravityModel::geometricGravityModel
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    gravityModel
    (
        typeName,
        dict,
        mesh
    ),
    ghf_
    (
        IOobject
        (
            "ghf",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("ghf_", dimAcceleration*dimLength, 0.0)
    ),
    hRef_
    (
        IOobject
        (
            "hRef",
            mesh.time().constant(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        dimensionedScalar(dimLength, Zero)
    )
{
    // Option to specify number of initial reconstruct calls.
    // Observed to be beneficial for plicRDF with high accuracy requirements.
    const label nInitReconstructs =
        dict.getOrDefault<label>("nInitReconstructs", 3);

    reconstructionSchemes& surf =
        mesh_.lookupObjectRef<reconstructionSchemes>("reconstructionScheme");

    for (label i = 0; i < nInitReconstructs; ++i)
    {
        surf.reconstruct(true);
    }

    calculateGhf();
}


// * * * * * * * * * * * * * * Public Access Member Functions  * * * * * * * //


void Foam::geometricGravityModel::calculateGhf()
{
    // Read only if mesh changed would be clever
    const uniformDimensionedVectorField& g =
        meshObjects::gravity::New(mesh_.time());

    dimensionedScalar ghRef
    (
        mag(g.value()) > SMALL
      ? g & (cmptMag(g.value())/mag(g.value()))*hRef_
      : dimensionedScalar("ghRef", g.dimensions()*dimLength, 0)
    );

    // Initial population of ghf - to be overwritten at fluid interface below
    ghf_ = (g & mesh_.Cf()) - ghRef;

    reconstructionSchemes& surf =
        mesh_.lookupObjectRef<reconstructionSchemes>("reconstructionScheme");

    surf.reconstruct(false);

    const boolList& isSurfCell = surf.interfaceCell();

    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
    scalarField& ghfIn = ghf_.primitiveFieldRef();

    // Setting ghf on internal faces of one or two interface cells
    forAll(ghfIn, facei)
    {
        if (isSurfCell[own[facei]])
        {
            point Ci = surf.centre()[own[facei]];
            if (isSurfCell[nei[facei]])
            {
                const scalar w = mesh_.weights()[facei];
                Ci = w*Ci + (1 - w)*surf.centre()[nei[facei]];
            }
            ghfIn[facei] = (g.value() & Ci) - ghRef.value();
        }
        else if (isSurfCell[nei[facei]])
        {
            const point Ci = surf.centre()[nei[facei]];
            ghfIn[facei] = (g.value() & Ci) - ghRef.value();
        }
    }

    // Setting ghf on boundary faces of interface cells
    forAll(patches, patchi)
    {
        scalarField& ghfbi = ghf_.boundaryFieldRef()[patchi];
        forAll(ghfbi, facei)
        {
            const label ownerCell = own[patches[patchi].start() + facei];
            if (isSurfCell[ownerCell])
                ghfbi[facei] =
                    (g.value() & surf.centre()[ownerCell]) - ghRef.value();
        }
    }


    // Declarations for interface-neighbour cell treatment further below.
    labelList nContrib(mesh_.nCells(), 0);
    boolList isNeighbourCell(mesh_.nCells(), false);
    scalarField ghNeiCell(mesh_.nCells(), 0.0);
    boolList faceRemoteIsInterface(mesh_.nBoundaryFaces(), false);

    // Single processor communication for interface-cell faces.
    // The receive loop simultaneously sets ghf on the face (interface-cell
    // pass) and accumulates the remote interface gh into ghNeiCell
    // (neighbour-cell pass), avoiding a second latency-dominated round over
    // the same data.
    if (Pstream::parRun())
    {
        DynamicList<label> neighProcs;
        PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

        for (const polyPatch& pp : patches)
        {
            const auto* ppp = isA<processorPolyPatch>(pp);
            if (ppp && pp.nPoints())
            {
                const auto& procPatch = *ppp;
                const label nbrProci = procPatch.neighbProcNo();
                neighProcs.append(nbrProci);
                UOPstream toNbr(nbrProci, pBufs);

                const scalarField& ghfbi =
                    ghf_.boundaryField()[procPatch.index()];
                List<label> interfaceFaces;
                List<scalar> interfaceFaceValues;
                forAll(ghfbi, facei)
                {
                    if (isSurfCell[own[procPatch.start() + facei]])
                    {
                        interfaceFaces.append(facei);
                        interfaceFaceValues.append(ghfbi[facei]);
                    }
                }
                toNbr << interfaceFaces << interfaceFaceValues;
            }
        }

        pBufs.finishedNeighbourSends(neighProcs);

        for (const polyPatch& pp : patches)
        {
            const auto* ppp = isA<processorPolyPatch>(pp);
            if (ppp && pp.nPoints())
            {
                const auto& procPatch = *ppp;
                const label nbrProci = procPatch.neighbProcNo();

                UIPstream fromNeighb(nbrProci, pBufs);
                List<label> neighbourInterfaceFaces;
                List<scalar> neighbourInterfaceFaceValues;
                fromNeighb >> neighbourInterfaceFaces >> neighbourInterfaceFaceValues;

                scalarField& ghfbi = ghf_.boundaryFieldRef()[procPatch.index()];
                const scalarField& procWeights =
                    mesh_.weights().boundaryField()[procPatch.index()];

                forAll(neighbourInterfaceFaces, fi)
                {
                    const label facei = neighbourInterfaceFaces[fi];
                    const label ownerCell = own[procPatch.start() + facei];
                    const scalar w = procWeights[facei];
                    if (isSurfCell[ownerCell])
                    {
                        // Both sides are interface cells: blend face value.
                        ghfbi[facei] =
                            w*ghfbi[facei] + (1 - w)*neighbourInterfaceFaceValues[fi];
                    }
                    else
                    {
                        // Remote is interface, local is not: set face value
                        // and accumulate remote interface gh into ghNeiCell.
                        ghfbi[facei] = neighbourInterfaceFaceValues[fi];
                        isNeighbourCell[ownerCell] = true;
                        ghNeiCell[ownerCell] += neighbourInterfaceFaceValues[fi];
                        nContrib[ownerCell]++;
                        const label bFacei =
                            procPatch.start() + facei - mesh_.nInternalFaces();
                        faceRemoteIsInterface[bFacei] = true;
                    }
                }
            }
        }
    }

    // Cyclic patches: same face-value logic as internal faces.
    forAll(patches, patchi)
    {
        const fvPatch& patch = mesh_.boundary()[patchi];
        if (!isA<cyclicFvPatch>(patch)) continue;

        const cyclicFvPatch& cyclicPatch =
            refCast<const cyclicFvPatch>(patch);
        if (cyclicPatch.neighbPatchID() < patchi) continue;

        scalarField& patchField = ghf_.boundaryFieldRef()[patchi];
        scalarField& neighbPatchField =
            ghf_.boundaryFieldRef()[cyclicPatch.neighbPatchID()];

        const label patchStart = patches[patchi].start();
        const label neighbStart =
            patches[cyclicPatch.neighbPatchID()].start();

        forAll(patchField, facei)
        {
            const label cellA = own[patchStart + facei];
            const label cellB = own[neighbStart + facei];

            scalar ghf = 0;
            if (isSurfCell[cellA])
            {
                point Ci = surf.centre()[cellA];
                if (isSurfCell[cellB])
                {
                    const scalar w =
                        mesh_.weights().boundaryField()[patchi][facei];
                    Ci = w*Ci + (1 - w)*surf.centre()[cellB];
                }
                ghf = (g.value() & Ci) - ghRef.value();
                patchField[facei] = ghf;
                neighbPatchField[facei] = ghf;
            }
            else if (isSurfCell[cellB])
            {
                const point Ci = surf.centre()[cellB];
                ghf = (g.value() & Ci) - ghRef.value();
                patchField[facei] = ghf;
                neighbPatchField[facei] = ghf;
            }
        }
    }

    // At this point all faces of interface cells are set correctly.

    // Now we treat faces of cells that are face neighbours to interface cells
    // The intended behaviour is as follows:
    // For each cellI (face) neighbouring one or more interface cells calculate
    // their average interface height and use that for ghf on all faces of cellI
    // that are not on an interface cell.
    // So the procedure is to go through all interface cells and visit all its
    // neighbour cells that is not an interface cell itself and add ghf to them
    // mark them as interface neighbours and bump up a counter for them saying
    // how many interface cells they neighbour.

    // Accumulate gh for interface-neighbour cells from LOCAL internal
    // face-neighbours. Remote processor-boundary contributions were already
    // accumulated in the communication round above.
    {
        const cellList& cells = mesh_.cells();
        forAll(cells, celli)
        {
            if (!isSurfCell[celli]) continue;

            const scalar ghCelli =
                (g.value() & surf.centre()[celli]) - ghRef.value();
            forAll(cells[celli], fi)
            {
                const label facei = cells[celli][fi];
                if (facei >= mesh_.nInternalFaces()) continue;
                const label nbrCelli =
                    (own[facei] == celli) ? nei[facei] : own[facei];
                if (!isSurfCell[nbrCelli])
                {
                    isNeighbourCell[nbrCelli] = true;
                    ghNeiCell[nbrCelli] += ghCelli;
                    nContrib[nbrCelli]++;
                }
            }
        }
    }


    // Normalise ghNeiCell (local + processor-boundary contributions combined)
    forAll(ghNeiCell, celli)
    {
        if (nContrib[celli] > 0)
            ghNeiCell[celli] = ghNeiCell[celli] / nContrib[celli];
    }

    // Now copy ghNeiCell out to all faces that do not have an interface cell
    // on the other side.
    // isNeighbourCell can be discarded - just use nContrib > 0.
    // ghNeiCell must also be copied to faces on processor patches if there is
    // not an interface cell on the other side.
    // This can be checked locally: If nContrib == nLocalSurroundingInterfaceCells
    // Then then there is NOT an interface on the other side of the processor
    // patch and we must copy ghNeiCell to the face patch.
    // A face can be shared between two neighbour-to-interface cells. In that
    // case we must use the weighted average value on the face. This can also
    // happen for a face on a processor patch. In that case we must pass values
    // across the patch again.


    // Setting ghf on internal faces of one or two interface neighbour cells
    forAll(ghfIn, facei)
    {
        if (isSurfCell[own[facei]] || isSurfCell[nei[facei]]) continue;
        if (isNeighbourCell[own[facei]])
        {
            ghfIn[facei] = ghNeiCell[own[facei]];
            if (isNeighbourCell[nei[facei]])
            {
                const scalar w = mesh_.weights()[facei];
                ghfIn[facei] = w*ghfIn[facei] + (1 - w)*ghNeiCell[nei[facei]];
            }
        }
        else if (isNeighbourCell[nei[facei]])
        {
            ghfIn[facei] = ghNeiCell[nei[facei]];
        }
    }

    // Setting ghf on boundary faces of interface neighbour cells
    forAll(patches, patchi)
    {
        if (isA<cyclicFvPatch>(mesh_.boundary()[patchi])) continue;
        if (isA<processorPolyPatch>(patches[patchi])) continue;

        scalarField& ghfbi = ghf_.boundaryFieldRef()[patchi];
        forAll(ghfbi, facei)
        {
            const label ownerCell = own[patches[patchi].start() + facei];
            if (isNeighbourCell[ownerCell])
                ghfbi[facei] = ghNeiCell[ownerCell];
        }
    }

    // Setting ghf on cyclic patch faces of interface neighbour cells.
    // Process each pair once (skip the patch with the higher index).
    forAll(patches, patchi)
    {
        const fvPatch& patch = mesh_.boundary()[patchi];
        if (!isA<cyclicFvPatch>(patch)) continue;

        const cyclicFvPatch& cyclicPatch =
            refCast<const cyclicFvPatch>(patch);
        if (cyclicPatch.neighbPatchID() < patchi) continue;

        scalarField& patchField = ghf_.boundaryFieldRef()[patchi];
        scalarField& neighbPatchField =
            ghf_.boundaryFieldRef()[cyclicPatch.neighbPatchID()];

        const label patchStart = patches[patchi].start();
        const label neighbStart =
            patches[cyclicPatch.neighbPatchID()].start();

        forAll(patchField, facei)
        {
            const label cellA = own[patchStart + facei];
            const label cellB = own[neighbStart + facei];

            if (isSurfCell[cellA] || isSurfCell[cellB]) continue;
            if (!isNeighbourCell[cellA] && !isNeighbourCell[cellB]) continue;

            scalar fv;
            if (isNeighbourCell[cellA] && isNeighbourCell[cellB])
            {
                const scalar w =
                    mesh_.weights().boundaryField()[patchi][facei];
                fv = w*ghNeiCell[cellA] + (1 - w)*ghNeiCell[cellB];
            }
            else if (isNeighbourCell[cellA])
                fv = ghNeiCell[cellA];
            else
                fv = ghNeiCell[cellB];

            patchField[facei] = fv;
            neighbPatchField[facei] = fv;
        }
    }


    // Setting ghf on processor patch faces of interface neighbour cells.
    // Faces where the remote was an interface cell are already correct from
    // pass 1 (faceRemoteIsInterface guards those). For all other neighbour-cell
    // faces the pre-step seeds ghfbi = ghNeiCell[local], which is the final
    // value when the remote is neither. The send/receive then overwrites with
    // the correctly blended value when both sides are neighbour cells,
    // matching the weighted-average treatment of ordinary internal faces.
    if (Pstream::parRun())
    {
        for (const polyPatch& pp : patches)
        {
            const auto* ppp = isA<processorPolyPatch>(pp);
            if (!ppp || !pp.nPoints()) continue;
            const auto& procPatch = *ppp;
            scalarField& ghfbi = ghf_.boundaryFieldRef()[procPatch.index()];
            forAll(ghfbi, facei)
            {
                const label ownerCell = own[procPatch.start() + facei];
                const label bFacei =
                    procPatch.start() + facei - mesh_.nInternalFaces();
                if
                (
                    isNeighbourCell[ownerCell]
                 && !isSurfCell[ownerCell]
                 && !faceRemoteIsInterface[bFacei]
                )
                    ghfbi[facei] = ghNeiCell[ownerCell];
            }
        }

        DynamicList<label> neighProcs;
        PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

        for (const polyPatch& pp : patches)
        {
            const auto* ppp = isA<processorPolyPatch>(pp);
            if (!ppp || !pp.nPoints()) continue;
            const auto& procPatch = *ppp;
            neighProcs.append(procPatch.neighbProcNo());
            UOPstream toNbr(procPatch.neighbProcNo(), pBufs);

            List<label> neiCellFaces;
            List<scalar> neiCellValues;
            forAll(ghf_.boundaryField()[procPatch.index()], facei)
            {
                const label ownerCell = own[procPatch.start() + facei];
                const label bFacei =
                    procPatch.start() + facei - mesh_.nInternalFaces();
                if
                (
                    isNeighbourCell[ownerCell]
                 && !isSurfCell[ownerCell]
                 && !faceRemoteIsInterface[bFacei]
                )
                {
                    neiCellFaces.append(facei);
                    neiCellValues.append(ghNeiCell[ownerCell]);
                }
            }
            toNbr << neiCellFaces << neiCellValues;
        }

        pBufs.finishedNeighbourSends(neighProcs);

        for (const polyPatch& pp : patches)
        {
            const auto* ppp = isA<processorPolyPatch>(pp);
            if (!ppp || !pp.nPoints()) continue;
            const auto& procPatch = *ppp;

            UIPstream fromNeighb(procPatch.neighbProcNo(), pBufs);
            List<label> remoteNeiCellFaces;
            List<scalar> remoteNeiCellValues;
            fromNeighb >> remoteNeiCellFaces >> remoteNeiCellValues;

            scalarField& ghfbi = ghf_.boundaryFieldRef()[procPatch.index()];
            const scalarField& procWeights =
                mesh_.weights().boundaryField()[procPatch.index()];

            forAll(remoteNeiCellFaces, fi)
            {
                const label facei = remoteNeiCellFaces[fi];
                const label ownerCell = own[procPatch.start() + facei];
                if (isSurfCell[ownerCell]) continue;
                const scalar w = procWeights[facei];
                if (isNeighbourCell[ownerCell])
                    ghfbi[facei] =
                        w*ghNeiCell[ownerCell] + (1 - w)*remoteNeiCellValues[fi];
                else
                    ghfbi[facei] = remoteNeiCellValues[fi];
            }
        }
    }
}


Foam::tmp<Foam::surfaceScalarField> Foam::geometricGravityModel::force()
{
    calculateGhf();
    const volScalarField& rho = mesh_.lookupObject<volScalarField>("rho");
    return -ghf_*fvc::snGrad(rho);
}

// ************************************************************************* //
