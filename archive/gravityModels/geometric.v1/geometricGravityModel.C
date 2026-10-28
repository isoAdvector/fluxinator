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

    ghf_ = (g & mesh_.Cf()) - ghRef;

    reconstructionSchemes& surf =
        mesh_.lookupObjectRef<reconstructionSchemes>("reconstructionScheme");

/*
    volVectorField nCopy = surf.normal();
    volVectorField cCopy = surf.centre();
*/
    // Force update of interface geometry data.
    // This has bee found to be necessary and make a real difference.
    // The solver should do this immediately after alpha update.
     surf.reconstruct(true);

//    surf.reconstruct(false);

/*
    volScalarField normalDiff = mag(surf.normal() - nCopy);
    volScalarField centreDiff = mag(surf.centre() - cCopy);
    Info << "After surf.reconstruct(true):" << endl;
    Info << "Max(|n-nold|) = " << max(normalDiff) << endl;
    Info << "Max(|c-cold|) = " << max(centreDiff) << endl;
*/

    const DynamicField<label>& interfaceLabels = surf.interfaceLabels();
    boolList isSurfCell(mesh_.nCells(), false);
    forAll(interfaceLabels, li)
    {
        isSurfCell[interfaceLabels[li]] = true;
    }

    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();

    // Setting internal ghf values
    scalarField& ghfIn = ghf_.primitiveFieldRef();
    forAll(ghfIn, facei)
    {
        const label ownerCell(own[facei]);
        const label neiCell(nei[facei]);
        point Ci(Zero);
        if (isSurfCell[ownerCell])
        {
            Ci = surf.centre()[ownerCell];
            if (isSurfCell[neiCell])
            {
                // Simple average
                // Ci += surf.centre()[neiCell];
                // Ci *= 0.5;
                // Distance-weighted interpolation (owner weight from mesh)
                const scalar w = mesh_.weights()[facei];
                Ci = w*surf.centre()[ownerCell] + (1 - w)*surf.centre()[neiCell];
            }
        }
        else if (isSurfCell[neiCell])
        {
            Ci = surf.centre()[neiCell];
        }
        else
        {
            Ci = mesh_.Cf()[facei];
        }
        ghfIn[facei] = (g.value() & Ci) -  ghRef.value();
    }


    // Setting boundary ghf values
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    forAll(patches, patchi)
    {
        scalarField& ghfbi = ghf_.boundaryFieldRef()[patchi];

        forAll(ghfbi, facei)
        {
            const label globalFacei = facei + patches[patchi].start();
            const label ownerCell(own[globalFacei]);

            point Ci(Zero);
            if (isSurfCell[ownerCell])
            {
                Ci = surf.centre()[ownerCell];
            }
            else
            {
                Ci = mesh_.Cf().boundaryField()[patchi][facei];
            }
            ghfbi[facei] = (g.value() & Ci) - ghRef.value();
        }
    }

    // Correcting ghf on cyclic patches
    // Loop over all patches
    forAll(patches, patchi)
    {
        const fvPatch& patch = mesh_.boundary()[patchi];

        // Check if the patch is cyclic
        if (isA<cyclicFvPatch>(patch))
        {
    //        Info << "Patch " << patchi << " is cyclic" << endl;
            const cyclicFvPatch& cyclicPatch =
                refCast<const cyclicFvPatch>(patch);

            // Access field values on the patches
            scalarField& patchField = ghf_.boundaryFieldRef()[patchi];
            scalarField& neighbPatchField =
                ghf_.boundaryFieldRef()[cyclicPatch.neighbPatchID()];

            // Set field values to the weighted mean of the corresponding face pairs
            forAll(patchField, facei)
            {
                // Simple average
                // scalar faceValue =
                //     0.5 * (patchField[facei] + neighbPatchField[facei]);
                // Distance-weighted interpolation (owner weight from mesh)
                const scalar w =
                    mesh_.weights().boundaryField()[patchi][facei];
                scalar faceValue =
                    w*patchField[facei] + (1 - w)*neighbPatchField[facei];
                patchField[facei] = faceValue;
                neighbPatchField[facei] = faceValue;
            }
        }
    }


    // For each non-interface cell, if it has one or more face-neighbours that
    // are interface cells, compute gh as the average of ghf_ on the shared
    // faces between this cell and its interface-cell face-neighbours.
    scalarField ghCell(mesh_.nCells(), 0.0);
    boolList isNeighbourCell(mesh_.nCells(), false);

    const cellList& cells = mesh_.cells();

    forAll(cells, celli)
    {
        if (isSurfCell[celli]) continue;

        scalar ghSum = 0.0;
        label nSharedFaces = 0;

        forAll(cells[celli], fi)
        {
            const label facei = cells[celli][fi];

            if (facei >= mesh_.nInternalFaces()) continue;

            const label nbrCelli =
                (own[facei] == celli) ? nei[facei] : own[facei];

            if (isSurfCell[nbrCelli])
            {
                isNeighbourCell[celli] = true;
                ghSum += ghfIn[facei];
                nSharedFaces++;
            }
        }

        if (isNeighbourCell[celli])
        {
            ghCell[celli] = ghSum / nSharedFaces;
        }
    }

    // Set ghf_ on internal faces that belong to interface-neighbour cells but
    // not to interface cells.
    forAll(ghfIn, facei)
    {
        const label ownerCell(own[facei]);
        const label neiCell(nei[facei]);

        if (!isSurfCell[ownerCell] && !isSurfCell[neiCell])
        {
            if (isNeighbourCell[ownerCell] && isNeighbourCell[neiCell])
            {
                const scalar w = mesh_.weights()[facei];
                ghfIn[facei] = w*ghCell[ownerCell] + (1 - w)*ghCell[neiCell];
            }
            else if (isNeighbourCell[ownerCell])
            {
                ghfIn[facei] = ghCell[ownerCell];
            }
            else if (isNeighbourCell[neiCell])
            {
                ghfIn[facei] = ghCell[neiCell];
            }
        }
    }

    // Set ghf_ on boundary faces of interface-neighbour cells.
    // Cyclic patches are corrected below; processor patches are set here so
    // the values are ready to be picked up by the processor send loop.
    forAll(patches, patchi)
    {
        if (isA<cyclicFvPatch>(mesh_.boundary()[patchi])) continue;

        scalarField& ghfbi = ghf_.boundaryFieldRef()[patchi];

        forAll(ghfbi, facei)
        {
            const label globalFacei = facei + patches[patchi].start();
            const label ownerCell(own[globalFacei]);

            if (!isSurfCell[ownerCell] && isNeighbourCell[ownerCell])
            {
                ghfbi[facei] = ghCell[ownerCell];
            }
        }
    }

    // Correct ghf_ on cyclic patches using the same decision logic as for
    // internal faces, covering interface cells and interface-neighbour cells.
    // This supersedes the earlier cyclic correction above.
    forAll(patches, patchi)
    {
        const fvPatch& patch = mesh_.boundary()[patchi];

        if (!isA<cyclicFvPatch>(patch)) continue;

        const cyclicFvPatch& cyclicPatch =
            refCast<const cyclicFvPatch>(patch);

        scalarField& patchField = ghf_.boundaryFieldRef()[patchi];
        scalarField& neighbPatchField =
            ghf_.boundaryFieldRef()[cyclicPatch.neighbPatchID()];

        const label patchStart = patches[patchi].start();
        const label neighbPatchStart =
            patches[cyclicPatch.neighbPatchID()].start();

        forAll(patchField, facei)
        {
            const label ownerCell = own[patchStart + facei];
            const label nbrOwnerCell = own[neighbPatchStart + facei];

            scalar faceValue;

            if (isSurfCell[ownerCell] && isSurfCell[nbrOwnerCell])
            {
                const scalar w =
                    mesh_.weights().boundaryField()[patchi][facei];
                faceValue =
                    (g.value() &
                        ( w*surf.centre()[ownerCell]
                        + (1 - w)*surf.centre()[nbrOwnerCell]))
                    - ghRef.value();
            }
            else if (isSurfCell[ownerCell])
            {
                faceValue =
                    (g.value() & surf.centre()[ownerCell]) - ghRef.value();
            }
            else if (isSurfCell[nbrOwnerCell])
            {
                faceValue =
                    (g.value() & surf.centre()[nbrOwnerCell]) - ghRef.value();
            }
            else if (isNeighbourCell[ownerCell] && isNeighbourCell[nbrOwnerCell])
            {
                const scalar w =
                    mesh_.weights().boundaryField()[patchi][facei];
                faceValue =
                    w*ghCell[ownerCell] + (1 - w)*ghCell[nbrOwnerCell];
            }
            else if (isNeighbourCell[ownerCell])
            {
                faceValue = ghCell[ownerCell];
            }
            else if (isNeighbourCell[nbrOwnerCell])
            {
                faceValue = ghCell[nbrOwnerCell];
            }
            else
            {
                continue;
            }

            patchField[facei] = faceValue;
            neighbPatchField[facei] = faceValue;
        }
    }


    // Setting ghf on processor patches
    if (Pstream::parRun())
    {
        DynamicList<label> neighProcs;
        PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

        // Send
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
                List<label> nextToInterfaceFaces;
                List<scalar> nextToInterfaceValues;
                forAll(ghfbi, facei)
                {
                    const label globalFacei = facei + procPatch.start();
                    const label ownerCell(own[globalFacei]);
                    if (isSurfCell[ownerCell])
                    {
                        interfaceFaces.append(facei);
                        interfaceFaceValues.append(ghfbi[facei]);
                    }
                    else if (isNeighbourCell[ownerCell])
                    {
                        nextToInterfaceFaces.append(facei);
                        nextToInterfaceValues.append(ghfbi[facei]);
                    }
                }

                toNbr << interfaceFaces << interfaceFaceValues
                       << nextToInterfaceFaces << nextToInterfaceValues;
            }
        }

        // Limited to involved neighbour procs
        pBufs.finishedNeighbourSends(neighProcs);


        // Receive and combine
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
                List<label> neighbourNextToInterfaceFaces;
                List<scalar> neighbourNextToInterfaceValues;

                fromNeighb >> neighbourInterfaceFaces
                    >> neighbourInterfaceFaceValues
                    >> neighbourNextToInterfaceFaces
                    >> neighbourNextToInterfaceValues;

                scalarField& ghfbi = ghf_.boundaryFieldRef()[procPatch.index()];

                const scalarField& procWeights =
                    mesh_.weights().boundaryField()[procPatch.index()];

                forAll(neighbourInterfaceFaces, fi)
                {
                    const label facei = neighbourInterfaceFaces[fi];
                    const label globalFacei = facei + procPatch.start();
                    const label ownerCell(own[globalFacei]);
                    if (isSurfCell[ownerCell])
                    {
                        // Simple average
                        // ghfbi[facei] += neighbourInterfaceFaceValues[fi];
                        // ghfbi[facei] *= 0.5;
                        // Distance-weighted interpolation (owner weight from mesh)
                        const scalar w = procWeights[facei];
                        ghfbi[facei] =
                            w*ghfbi[facei]
                          + (1 - w)*neighbourInterfaceFaceValues[fi];
                    }
                    else
                    {
                        ghfbi[facei] = neighbourInterfaceFaceValues[fi];
                    }
                }

                forAll(neighbourNextToInterfaceFaces, fi)
                {
                    const label facei = neighbourNextToInterfaceFaces[fi];
                    const label globalFacei = facei + procPatch.start();
                    const label ownerCell(own[globalFacei]);
                    if (isSurfCell[ownerCell]) continue;
                    const scalar w = procWeights[facei];
                    if (isNeighbourCell[ownerCell])
                    {
                        ghfbi[facei] =
                            w*ghfbi[facei]
                          + (1 - w)*neighbourNextToInterfaceValues[fi];
                    }
                    else
                    {
                        ghfbi[facei] = neighbourNextToInterfaceValues[fi];
                    }
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
