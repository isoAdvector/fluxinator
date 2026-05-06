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
    const boolList& isSurfCell = surf.interfaceCell();

    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
    scalarField& ghfIn = ghf_.primitiveFieldRef();

    // Accumulate gh for interface-neighbour cells from LOCAL face-neighbours.
    // Processor-boundary contributions are added during pass 1 receive below.
    // Raw sums are kept; normalisation happens after both contributions.
    scalarField ghSum(mesh_.nCells(), 0.0);
    labelList nContrib(mesh_.nCells(), 0);
    boolList isNeighbourCell(mesh_.nCells(), false);
    scalarField ghCell(mesh_.nCells(), 0.0);
    {
        const cellList& cells = mesh_.cells();
        forAll(cells, celli)
        {
            if (isSurfCell[celli]) continue;
            forAll(cells[celli], fi)
            {
                const label facei = cells[celli][fi];
                if (facei >= mesh_.nInternalFaces()) continue;
                const label nbrCelli =
                    (own[facei] == celli) ? nei[facei] : own[facei];
                if (isSurfCell[nbrCelli])
                {
                    isNeighbourCell[celli] = true;
                    ghSum[celli] +=
                        (g.value() & surf.centre()[nbrCelli]) - ghRef.value();
                    nContrib[celli]++;
                }
            }
        }
    }

    // Set processor boundary face values for interface cells so they are
    // ready to send in pass 1.
    if (Pstream::parRun())
    {
        for (const polyPatch& pp : patches)
        {
            const auto* ppp = isA<processorPolyPatch>(pp);
            if (ppp && pp.nPoints())
            {
                const auto& procPatch = *ppp;
                scalarField& ghfbi = ghf_.boundaryFieldRef()[procPatch.index()];
                forAll(ghfbi, facei)
                {
                    const label ownerCell = own[procPatch.start() + facei];
                    if (isSurfCell[ownerCell])
                        ghfbi[facei] =
                            (g.value() & surf.centre()[ownerCell])
                          - ghRef.value();
                }
            }
        }
    }

    // Processor pass 1: communicate interface cell face values.
    // The receive loop also marks isNeighbourCell and accumulates ghCell for
    // local cells that are face-adjacent to remote interface cells.
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
                        ghfbi[facei] =
                            w*ghfbi[facei] + (1 - w)*neighbourInterfaceFaceValues[fi];
                    }
                    else
                    {
                        ghfbi[facei] = neighbourInterfaceFaceValues[fi];
                        isNeighbourCell[ownerCell] = true;
                        ghSum[ownerCell] += ghfbi[facei];
                        nContrib[ownerCell]++;
                    }
                }
            }
        }
    }

    // Normalise ghCell (local + processor-boundary contributions combined)
    forAll(ghCell, celli)
    {
        if (nContrib[celli] > 0)
            ghCell[celli] = ghSum[celli] / nContrib[celli];
    }

    // Returns the ghf value for a face between cellA and cellB with owner
    // weight w. Priority: interface cell > interface-neighbour cell > fallback.
    auto faceGhf =
        [&](label cellA, label cellB, scalar w, scalar fallback) -> scalar
    {
        if (isSurfCell[cellA] || isSurfCell[cellB])
        {
            const point Ci =
                (isSurfCell[cellA] && isSurfCell[cellB])
              ? w*surf.centre()[cellA] + (1 - w)*surf.centre()[cellB]
              : isSurfCell[cellA] ? surf.centre()[cellA] : surf.centre()[cellB];
            return (g.value() & Ci) - ghRef.value();
        }
        if (isNeighbourCell[cellA] && isNeighbourCell[cellB])
            return w*ghCell[cellA] + (1 - w)*ghCell[cellB];
        if (isNeighbourCell[cellA]) return ghCell[cellA];
        if (isNeighbourCell[cellB]) return ghCell[cellB];
        return fallback;
    };

    // Internal faces
    forAll(ghfIn, facei)
    {
        ghfIn[facei] = faceGhf
        (
            own[facei], nei[facei],
            mesh_.weights()[facei],
            (g.value() & mesh_.Cf()[facei]) - ghRef.value()
        );
    }

    // Non-cyclic, non-processor boundary faces
    forAll(patches, patchi)
    {
        if (isA<cyclicFvPatch>(mesh_.boundary()[patchi])) continue;
        if (isA<processorPolyPatch>(patches[patchi])) continue;

        scalarField& ghfbi = ghf_.boundaryFieldRef()[patchi];
        const vectorField& Cfb = mesh_.Cf().boundaryField()[patchi];

        forAll(ghfbi, facei)
        {
            const label ownerCell = own[patches[patchi].start() + facei];
            if (isSurfCell[ownerCell])
                ghfbi[facei] =
                    (g.value() & surf.centre()[ownerCell]) - ghRef.value();
            else if (isNeighbourCell[ownerCell])
                ghfbi[facei] = ghCell[ownerCell];
            else
                ghfbi[facei] = (g.value() & Cfb[facei]) - ghRef.value();
        }
    }

    // Cyclic patches: same face-value logic as internal faces.
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

            if ( !isSurfCell[cellA] && !isSurfCell[cellB]
              && !isNeighbourCell[cellA] && !isNeighbourCell[cellB]) continue;

            const scalar w =
                mesh_.weights().boundaryField()[patchi][facei];
            const scalar fv = faceGhf(cellA, cellB, w, 0);

            patchField[facei] = fv;
            neighbPatchField[facei] = fv;
        }
    }

    // Set processor boundary face values for neighbour cells so they are
    // ready to send in pass 2.
    if (Pstream::parRun())
    {
        for (const polyPatch& pp : patches)
        {
            const auto* ppp = isA<processorPolyPatch>(pp);
            if (ppp && pp.nPoints())
            {
                const auto& procPatch = *ppp;
                scalarField& ghfbi = ghf_.boundaryFieldRef()[procPatch.index()];
                forAll(ghfbi, facei)
                {
                    const label ownerCell = own[procPatch.start() + facei];
                    if (!isSurfCell[ownerCell] && isNeighbourCell[ownerCell])
                        ghfbi[facei] = ghCell[ownerCell];
                }
            }
        }
    }

    // Processor pass 2: communicate neighbour cell face values.
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
                List<label> nextToInterfaceFaces;
                List<scalar> nextToInterfaceValues;
                forAll(ghfbi, facei)
                {
                    const label ownerCell = own[procPatch.start() + facei];
                    if (!isSurfCell[ownerCell] && isNeighbourCell[ownerCell])
                    {
                        nextToInterfaceFaces.append(facei);
                        nextToInterfaceValues.append(ghfbi[facei]);
                    }
                }
                toNbr << nextToInterfaceFaces << nextToInterfaceValues;
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
                List<label> neighbourNextToInterfaceFaces;
                List<scalar> neighbourNextToInterfaceValues;
                fromNeighb >> neighbourNextToInterfaceFaces
                    >> neighbourNextToInterfaceValues;

                scalarField& ghfbi = ghf_.boundaryFieldRef()[procPatch.index()];
                const scalarField& procWeights =
                    mesh_.weights().boundaryField()[procPatch.index()];

                forAll(neighbourNextToInterfaceFaces, fi)
                {
                    const label facei = neighbourNextToInterfaceFaces[fi];
                    const label ownerCell = own[procPatch.start() + facei];
                    if (isSurfCell[ownerCell]) continue;
                    const scalar w = procWeights[facei];
                    if (isNeighbourCell[ownerCell])
                        ghfbi[facei] =
                            w*ghfbi[facei] + (1 - w)*neighbourNextToInterfaceValues[fi];
                    else
                        ghfbi[facei] = neighbourNextToInterfaceValues[fi];
                }
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
