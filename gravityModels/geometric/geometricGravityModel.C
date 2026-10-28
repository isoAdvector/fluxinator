/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2025 Henning Scheufler
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
#include "syncTools.H"


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
    gh_
    (
        IOobject
        (
            "gh",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("gh_", dimAcceleration*dimLength, 0.0)
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
    ),
    nLayers_(dict.getOrDefault<label>("nLayers", 2))
{
    if (nLayers_ < 0)
    {
        FatalIOErrorInFunction(dict)
            << "nLayers = " << nLayers_ << " must be zero or positive"
            << exit(FatalIOError);
    }

    Info<< "geometric gravity model: nLayers = " << nLayers_ << endl;

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

    const fvMesh& mesh = gh_.mesh();
    zoneDistribute& exchangeFields = zoneDistribute::New(mesh);

    const uniformDimensionedVectorField& g =
        meshObjects::gravity::New(mesh_.time());

    dimensionedScalar ghRef
    (
        mag(g.value()) > SMALL
      ? g & (cmptMag(g.value())/mag(g.value()))*hRef_
      : dimensionedScalar("ghRef", g.dimensions()*dimLength, 0)
    );

    gh_ = (g & mesh.C()) - ghRef;
    ghf_ = ((g & mesh.Cf()) - ghRef);

    if(mesh.foundObject<reconstructionSchemes>("reconstructionScheme"))
    {
        reconstructionSchemes& surf =
            mesh.lookupObjectRef<reconstructionSchemes>("reconstructionScheme");

//        surf.reconstruct();
        // Interface centres, normals and cells of the reconstruction, extended
        // below by faces that are themselves the interface
        volVectorField faceCentre("ghInterfaceCentre", surf.centre());
        volVectorField faceNormal("ghInterfaceNormal", surf.normal());
        boolList interfaceCell(surf.interfaceCell());

        // An interface lying exactly on mesh faces (e.g. a still water level
        // initialised on a face) leaves volume fractions of exactly 0 and 1
        // and no interface cells. Such a face separates two non-interface
        // cells with densities on either side of the mean. Both cells get the
        // face as their interface (area-weighted over several such faces).
        if (mesh.foundObject<volScalarField>("rho"))
        {
            const volScalarField& rho = mesh.lookupObject<volScalarField>("rho");
            const scalar rhoMid = 0.5*(gMin(rho.primitiveField()) + gMax(rho.primitiveField()));

            const labelUList& own = mesh.owner();
            const labelUList& nei = mesh.neighbour();
            const label nBFaces = mesh.nBoundaryFaces();

            boolList nbrInterfaceCell;
            syncTools::swapBoundaryCellList(mesh, interfaceCell, nbrInterfaceCell);
            scalarField nbrRho;
            syncTools::swapBoundaryCellList(mesh, rho.primitiveField(), nbrRho);

            vectorField sumCf(mesh.nCells(), Zero);
            vectorField sumSf(mesh.nCells(), Zero);
            scalarField sumA(mesh.nCells(), Zero);

            // Adds face facei to cell celli; Sf is oriented towards the lighter side
            const auto addFace = [&](label celli, label facei, scalar sign)
            {
                const scalar A = mesh.magSf()[facei];
                sumCf[celli] += A*mesh.Cf()[facei];
                sumSf[celli] += sign*mesh.Sf()[facei];
                sumA[celli] += A;
            };

            forAll(nei, facei)
            {
                const label P = own[facei], N = nei[facei];
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
                if (!pp.coupled()) continue;

                forAll(pp, i)
                {
                    const label facei = pp.start() + i;
                    const label bFacei = facei - mesh.nInternalFaces();
                    const label P = own[facei];
                    if
                    (
                        bFacei < nBFaces
                     && !interfaceCell[P] && !nbrInterfaceCell[bFacei]
                     && (rho[P] - rhoMid)*(nbrRho[bFacei] - rhoMid) < 0
                    )
                    {
                        const scalar sign = rho[P] > nbrRho[bFacei] ? 1 : -1;
                        const scalar A = mag(pp.faceAreas()[i]);
                        sumCf[P] += A*pp.faceCentres()[i];
                        sumSf[P] += sign*pp.faceAreas()[i];
                        sumA[P] += A;
                    }
                }
            }

            forAll(sumA, celli)
            {
                if (sumA[celli] > 0 && mag(sumSf[celli]) > 0)
                {
                    faceCentre[celli] = sumCf[celli]/sumA[celli];
                    faceNormal[celli] = sumSf[celli];
                    interfaceCell[celli] = true;
                }
            }
        }

        // Interface cells and nLayers rings of point neighbours, with the
        // ring number of each cell (0 for interface cells, -1 outside)
        boolList nearInterface(mesh.nCells(), false);
        volScalarField cellDistLevel
        (
            IOobject
            (
                "cellDistLevel",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            mesh,
            dimensionedScalar(dimless, -1)
        );

        markInterfaceRegion markIF(mesh);
        markIF.markCellsNearSurf
        (
            interfaceCell,
            nLayers_,
            nearInterface,
            cellDistLevel
        );

        const auto inLayer = [&cellDistLevel](label celli, label layer)
        {
            return mag(cellDistLevel[celli] - layer) < 0.5;
        };

        exchangeFields.setUpCommforZone(nearInterface,true);

        Map<vector > mapCentres =
            exchangeFields.getDatafromOtherProc(nearInterface,faceCentre);
        Map<vector > mapNormal =
            exchangeFields.getDatafromOtherProc(nearInterface,faceNormal);

        const labelListList& stencil = exchangeFields.getStencil();

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

        // Interface cells and the first ring
        forAll(interfaceCell,celli)
        {
            if(mag(faceNormal[celli]) != 0)
            {
                gh_[celli] = (faceCentre[celli] & g.value()) - ghRef.value();
                assigned[celli] = 1;
            }
            else if(inLayer(celli, 1))
            {
                // the to the interface
                vector averageSurfP = vector::zero;
                scalar avgWeight = 0;
                const point p = mesh.C()[celli];

                forAll(stencil[celli],i)
                {
                    const label& gblIdx = stencil[celli][i];
                    vector n = -exchangeFields.getValue
                        (
                            faceNormal,
                            mapNormal,
                            gblIdx
                        );

                    if (mag(n) != 0)
                    {
                        n /= mag(n);
                        vector c = exchangeFields.getValue
                        (
                            faceCentre,
                            mapCentres,
                            gblIdx
                        );
                        vector distanceToIntSeg = (c - p);
                        scalar weight = 0;

                        if (mag(distanceToIntSeg) != 0)
                        {
                            distanceToIntSeg /= mag(distanceToIntSeg);
                            weight = sqr(mag(distanceToIntSeg & n));
                        }
                        else // exactly on the center
                        {
                            weight = 1;
                        }
                        averageSurfP += c * weight;
                        avgWeight += weight;
                    }
                }

                if (avgWeight != 0)
                {
                    averageSurfP /= avgWeight;
                    gh_[celli] = (averageSurfP & g.value()) - ghRef.value();
                    assigned[celli] = 1;
                }
            }
            else
            {
                // do nothing
            }

        }

        // Further rings: average of the values already assigned to the
        // neighbours. Each ring is one pass with its own parallel exchange,
        // and the new values only take effect after the pass.
        for (label layer = 2; layer <= nLayers_; ++layer)
        {
            Map<scalar> mapGh =
                exchangeFields.getDatafromOtherProc(nearInterface, gh_);
            Map<scalar> mapAssigned =
                exchangeFields.getDatafromOtherProc(nearInterface, assigned);

            scalarField ghNew(gh_.primitiveField());
            boolList newlyAssigned(mesh.nCells(), false);

            forAll(nearInterface, celli)
            {
                if (!inLayer(celli, layer))
                {
                    continue;
                }

                scalar ghSum = 0;
                label nAssigned = 0;

                forAll(stencil[celli], i)
                {
                    const label& gblIdx = stencil[celli][i];

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

        gh_.correctBoundaryConditions();
        ghf_ = fvc::interpolate(gh_);

    }
    else
    {
        Info << "notFound" << endl;
    }

}


Foam::tmp<Foam::surfaceScalarField> Foam::geometricGravityModel::force()
{
    calculateGhf();
    const volScalarField& rho = mesh_.lookupObject<volScalarField>("rho");
    return -ghf_*fvc::snGrad(rho);
}


void Foam::geometricGravityModel::setNLayers(const label n)
{
    if (n < 0)
    {
        FatalErrorInFunction
            << "nLayers = " << n << " must be zero or positive"
            << exit(FatalError);
    }

    nLayers_ = n;

    Info<< "geometric gravity model: nLayers = " << nLayers_ << endl;
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
