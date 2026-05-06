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
                forAll(ghfbi, facei)
                {
                    const label globalFacei = facei + procPatch.start();
                    const label ownerCell(own[globalFacei]);
                    if (isSurfCell[ownerCell])
                    {
                        interfaceFaces.append(facei);
                        interfaceFaceValues.append(ghfbi[facei]);
                    }
                }

                /*
                Pout << "Processor " << Pstream::myProcNo() <<
                    " sending interfaceFaces = " << interfaceFaces <<
                    " and interfaceFaceValues = " << interfaceFaceValues <<
                    " to neighbour processor = " << nbrProci << endl;
                */

                toNbr << interfaceFaces << interfaceFaceValues;
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

                fromNeighb >> neighbourInterfaceFaces
                    >> neighbourInterfaceFaceValues;

                /*
                Pout << "Processor " << Pstream::myProcNo() <<
                    " received  neighbourInterfaceFaces = " <<  neighbourInterfaceFaces <<
                    " and neighbourInterfaceFaceValues = " << neighbourInterfaceFaceValues <<
                    " from neighbour processor: " << nbrProci << endl;
                */

                scalarField& ghfbi = ghf_.boundaryFieldRef()[procPatch.index()];

                const scalarField& procWeights =
                    mesh_.weights().boundaryField()[procPatch.index()];

                /*
                Pout << "Processor patch " << procPatch.index()
                     << " weights range: ["
                     << min(procWeights) << ", " << max(procWeights) << "]"
                     << endl;
                */

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
