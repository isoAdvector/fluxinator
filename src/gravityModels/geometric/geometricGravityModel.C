/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2025 Henning Scheufler
    Copyright (C) 2025-2026 Johan Roenby
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
#include "markInterfaceRegion.H"
#include "syncTools.H"

using namespace Foam::fluxinatorVoF;


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(geometricGravityModel, 0);
    addToRunTimeSelectionTable(gravityModel, geometricGravityModel, components);
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::geometricGravityModel::addInterfaceFaces
(
    boolList& interfaceCell,
    volVectorField& interfaceCentre,
    volVectorField& interfaceNormal
) const
{
    // A face separates two non-interface cells with densities on either side
    // of the mean. Both cells get the face as their interface segment
    // (area-weighted over several such faces), with the normal towards the
    // lighter side.
    const fvMesh& mesh = mesh_;
    const volScalarField& rho = mesh.lookupObject<volScalarField>("rho");
    const scalar rhoMid =
        0.5*(gMin(rho.primitiveField()) + gMax(rho.primitiveField()));

    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();

    boolList nbrInterfaceCell;
    syncTools::swapBoundaryCellList(mesh, interfaceCell, nbrInterfaceCell);
    scalarField nbrRho;
    syncTools::swapBoundaryCellList(mesh, rho.primitiveField(), nbrRho);

    vectorField sumCf(mesh.nCells(), Zero);
    vectorField sumSf(mesh.nCells(), Zero);
    scalarField sumA(mesh.nCells(), Zero);

    // Adds face facei to cell celli, with Sf towards the lighter side
    const auto addFace = [&](label celli, label facei, scalar sign)
    {
        const scalar A = mag(mesh.faceAreas()[facei]);
        sumCf[celli] += A*mesh.faceCentres()[facei];
        sumSf[celli] += sign*mesh.faceAreas()[facei];
        sumA[celli] += A;
    };

    forAll(nei, facei)
    {
        const label P = own[facei];
        const label N = nei[facei];
        if
        (
            !interfaceCell[P] && !interfaceCell[N]
         && (rho[P] - rhoMid)*(rho[N] - rhoMid) < 0
        )
        {
            const scalar sign = rho[P] > rho[N] ? 1 : -1;
            addFace(P, facei, sign);
            addFace(N, facei, sign);
        }
    }

    for (const polyPatch& pp : mesh.boundaryMesh())
    {
        if (!pp.coupled())
        {
            continue;
        }

        forAll(pp, i)
        {
            const label facei = pp.start() + i;
            const label bFacei = facei - mesh.nInternalFaces();
            const label P = own[facei];
            if
            (
                !interfaceCell[P] && !nbrInterfaceCell[bFacei]
             && (rho[P] - rhoMid)*(nbrRho[bFacei] - rhoMid) < 0
            )
            {
                addFace(P, facei, rho[P] > nbrRho[bFacei] ? 1 : -1);
            }
        }
    }

    forAll(sumA, celli)
    {
        if (sumA[celli] > 0 && mag(sumSf[celli]) > 0)
        {
            interfaceCentre[celli] = sumCf[celli]/sumA[celli];
            interfaceNormal[celli] = sumSf[celli];
            interfaceCell[celli] = true;
        }
    }
}


void Foam::geometricGravityModel::setInterfaceCellsAndFirstLayer
(
    const volVectorField& interfaceCentre,
    const volVectorField& interfaceNormal,
    const labelList& cellLayer,
    const boolList& nearInterface,
    zoneDistribute& exchangeFields,
    const scalar ghRef,
    volScalarField& assigned
)
{
    const fvMesh& mesh = mesh_;
    const vector g = meshObjects::gravity::New(mesh.time()).value();

    const Map<vector> mapCentres =
        exchangeFields.getDatafromOtherProc(nearInterface, interfaceCentre);
    const Map<vector> mapNormals =
        exchangeFields.getDatafromOtherProc(nearInterface, interfaceNormal);

    const labelListList& stencil = exchangeFields.getStencil();

    forAll(cellLayer, celli)
    {
        if (mag(interfaceNormal[celli]) != 0)
        {
            // Interface cell: the height of its interface centroid
            gh_[celli] = (interfaceCentre[celli] & g) - ghRef;
            assigned[celli] = 1;
        }
        else if (cellLayer[celli] == 1)
        {
            // Layer 1: the mean of the neighbouring interface centroids,
            // weighted by (d & n)^2, with d the unit vector from the cell
            // centre to the centroid and n the interface normal
            vector centroidSum = Zero;
            scalar weightSum = 0;
            const point& C = mesh.C()[celli];

            for (const label gblIdx : stencil[celli])
            {
                vector n =
                    -exchangeFields.getValue(interfaceNormal, mapNormals, gblIdx);

                if (mag(n) != 0)
                {
                    n /= mag(n);
                    const vector c =
                        exchangeFields.getValue
                        (
                            interfaceCentre,
                            mapCentres,
                            gblIdx
                        );
                    vector d = c - C;
                    scalar weight = 0;

                    if (mag(d) != 0)
                    {
                        d /= mag(d);
                        weight = sqr(mag(d & n));
                    }
                    else
                    {
                        // The centroid is at the cell centre
                        weight = 1;
                    }
                    centroidSum += c*weight;
                    weightSum += weight;
                }
            }

            if (weightSum != 0)
            {
                centroidSum /= weightSum;
                gh_[celli] = (centroidSum & g) - ghRef;
                assigned[celli] = 1;
            }
        }
    }
}


void Foam::geometricGravityModel::setFurtherLayers
(
    const labelList& cellLayer,
    const boolList& nearInterface,
    zoneDistribute& exchangeFields,
    volScalarField& assigned
)
{
    const labelListList& stencil = exchangeFields.getStencil();

    // One pass per layer, with its own parallel exchange; the new values only
    // take effect after the pass
    for (label layer = 2; layer <= nLayers_; ++layer)
    {
        const Map<scalar> mapGh =
            exchangeFields.getDatafromOtherProc(nearInterface, gh_);
        const Map<scalar> mapAssigned =
            exchangeFields.getDatafromOtherProc(nearInterface, assigned);

        scalarField ghNew(gh_.primitiveField());
        boolList newlyAssigned(cellLayer.size(), false);

        forAll(cellLayer, celli)
        {
            if (cellLayer[celli] != layer)
            {
                continue;
            }

            scalar ghSum = 0;
            label nAssigned = 0;

            for (const label gblIdx : stencil[celli])
            {
                if
                (
                    exchangeFields.getValue(assigned, mapAssigned, gblIdx)
                  > 0.5
                )
                {
                    ghSum += exchangeFields.getValue(gh_, mapGh, gblIdx);
                    ++nAssigned;
                }
            }

            if (nAssigned > 0)
            {
                ghNew[celli] = ghSum/nAssigned;
                newlyAssigned[celli] = true;
            }
        }

        gh_.primitiveFieldRef() = ghNew;

        forAll(newlyAssigned, celli)
        {
            if (newlyAssigned[celli])
            {
                assigned[celli] = 1;
            }
        }
    }
}


// Sets g*h from the reconstructed interface: (1) take the interface centroids
// and normals of the reconstruction, adding faces that are themselves the
// interface, (2) mark the interface cells and the layers around them, (3) set
// g*h in the interface cells and the first layer, (4) then in each further
// layer, (5) interpolate to the faces. Elsewhere g*h is that of the cell
// centre and the face centre.
void Foam::geometricGravityModel::calculateGhf()
{
    const fvMesh& mesh = mesh_;
    const uniformDimensionedVectorField& g =
        meshObjects::gravity::New(mesh.time());
    const dimensionedScalar ghRef(this->ghRef());

    gh_ = (g & mesh.C()) - ghRef;
    ghf_ = ((g & mesh.Cf()) - ghRef);

    const auto* surfPtr =
        mesh.findObject<reconstructionSchemes>("reconstructionScheme");
    if (!surfPtr)
    {
        FatalErrorInFunction
            << "The geometric gravity model needs the interface reconstruction"
            << " of a geometric VoF solver (isoAdvector)" << exit(FatalError);
    }

    // (1)
    volVectorField interfaceCentre("ghInterfaceCentre", surfPtr->centre());
    volVectorField interfaceNormal("ghInterfaceNormal", surfPtr->normal());
    boolList interfaceCell(surfPtr->interfaceCell());
    addInterfaceFaces(interfaceCell, interfaceCentre, interfaceNormal);

    // (2)
    boolList nearInterface;
    labelList cellLayer;
    markInterfaceRegion(mesh).markCellsNearSurf
    (
        interfaceCell,
        nLayers_,
        nearInterface,
        cellLayer
    );

    zoneDistribute& exchangeFields = zoneDistribute::New(mesh);
    exchangeFields.setUpCommforZone(nearInterface, true);

    // 1 where gh_ holds an interface-based value
    volScalarField assigned
    (
        IOobject
        (
            "ghAssigned",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh,
        dimensionedScalar(dimless, Zero)
    );

    // (3), (4)
    setInterfaceCellsAndFirstLayer
    (
        interfaceCentre,
        interfaceNormal,
        cellLayer,
        nearInterface,
        exchangeFields,
        ghRef.value(),
        assigned
    );
    setFurtherLayers(cellLayer, nearInterface, exchangeFields, assigned);

    // (5)
    gh_.correctBoundaryConditions();
    ghf_ = fvc::interpolate(gh_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::geometricGravityModel::geometricGravityModel
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    gravityModel(typeName, dict, mesh),
    ghf_
    (
        IOobject
        (
            "ghf",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh,
        dimensionedScalar(dimAcceleration*dimLength, Zero)
    ),
    gh_
    (
        IOobject
        (
            "gh",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh,
        dimensionedScalar(dimAcceleration*dimLength, Zero)
    ),
    nLayers_(dict.getOrDefault<label>("nLayers", 2))
{
    if (nLayers_ < 0)
    {
        FatalIOErrorInFunction(dict)
            << "nLayers = " << nLayers_ << " must be zero or positive"
            << exit(FatalIOError);
    }

    Info<< "geometric gravity model: nLayers " << nLayers_ << endl;

    calculateGhf();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::surfaceScalarField> Foam::geometricGravityModel::force()
{
    calculateGhf();
    const volScalarField& rho = mesh_.lookupObject<volScalarField>("rho");
    return -ghf_*fvc::snGrad(rho);
}


Foam::tmp<Foam::surfaceScalarField> Foam::geometricGravityModel::ghf() const
{
    return ghf_;
}


Foam::tmp<Foam::volScalarField> Foam::geometricGravityModel::gh() const
{
    return gh_;
}


// ************************************************************************* //
