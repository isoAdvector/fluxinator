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

#include "densityWeightedGravityModel.H"
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
    defineTypeNameAndDebug(densityWeightedGravityModel, 0);
    addToRunTimeSelectionTable(gravityModel, densityWeightedGravityModel, components);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::densityWeightedGravityModel::densityWeightedGravityModel
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
    ghAssigned_
    (
        IOobject
        (
            "ghAssigned",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar(dimless, Zero)
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
    nLayers_(dict.getOrDefault<label>("nLayers", 1)),
    ghAtCellCentre_(dict.getOrDefault<bool>("ghAtCellCentre", false)),
    densityWeighted_(dict.getOrDefault<bool>("densityWeighted", false)),
    drhoTol_(dict.getOrDefault<scalar>("drhoTol", 1e-3))
{
    if (nLayers_ < 0)
    {
        FatalIOErrorInFunction(dict)
            << "nLayers = " << nLayers_ << " must be zero or positive"
            << exit(FatalIOError);
    }

    Info<< "densityWeighted gravity model: nLayers = " << nLayers_
        << ", ghAtCellCentre = " << Switch(ghAtCellCentre_)
        << ", densityWeighted = " << Switch(densityWeighted_) << endl;

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

void Foam::densityWeightedGravityModel::calculateGhf()
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
        const volVectorField& faceCentre = surf.centre();
        const volVectorField& faceNormal = surf.normal();

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
            surf.interfaceCell(),
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
        const vector gHat =
            mag(g.value()) > SMALL ? g.value()/mag(g.value()) : vector::zero;

        forAll(surf.interfaceCell(),celli)
        {
            if(mag(faceNormal[celli]) != 0)
            {
                // Height of the interface segment centroid, or of the point
                // where the reconstructed plane crosses the vertical through
                // the cell centre
                point interfaceP = faceCentre[celli];

                if (ghAtCellCentre_)
                {
                    const vector nHat = normalised(faceNormal[celli]);
                    const scalar nDotG = nHat & gHat;

                    // Skip where the interface is nearly parallel to gravity:
                    // the vertical through the cell centre then barely meets
                    // the plane and the projection is ill conditioned
                    if (mag(nDotG) > 0.01)
                    {
                        const point& Ci = mesh.C()[celli];
                        interfaceP =
                            Ci + (((faceCentre[celli] - Ci) & nHat)/nDotG)*gHat;
                    }
                }

                gh_[celli] = (interfaceP & g.value()) - ghRef.value();
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
        ghAssigned_ = assigned;
        ghAssigned_.correctBoundaryConditions();

        ghf_ = fvc::interpolate(gh_);

    }
    else
    {
        ghAssigned_ = dimensionedScalar(dimless, Zero);
        Info << "notFound" << endl;
    }

}


void Foam::densityWeightedGravityModel::densityWeightedForce
(
    const volScalarField& rho,
    surfaceScalarField& F
) const
{
    // On the faces that can carry it, take the force to be -snGrad(rho*gh)
    // instead of -ghf*snGrad(rho).
    //
    // The two agree wherever g*h is constant across the face, so a level
    // interface is unaffected. Where the interface slopes they differ, and
    // -snGrad(rho*gh) is the better choice: being the face gradient of a cell
    // field it is a discrete gradient, so its loop imbalance vanishes
    // identically and the pressure can balance it whatever the interface does
    // inside the cells.
    //
    // A corrected snGrad adds a non-orthogonal correction built from the Gauss
    // cell gradient, which samples rho*gh on every face of the cell and so
    // reaches the cell's face neighbours and, next to a physical boundary, the
    // boundary value of g*h - the height of the face centre rather than of the
    // interface. Substituting the force on a face whose gradient stencil
    // leaves the interface region mixes the two definitions of g*h and breaks
    // the balance. So a cell only qualifies if it and all of its face
    // neighbours carry an interface-based g*h, and a cell against a physical
    // boundary never qualifies. With an uncorrected snGrad the restriction
    // costs nothing that matters, and with a corrected one it is what keeps
    // the hydrostatic balance of Section 5 exact.
    const fvMesh& mesh = gh_.mesh();
    const label nInt = mesh.nInternalFaces();

    const surfaceScalarField snRhoGh(fvc::snGrad(rho*gh_));

    // Confine the substitution to the faces that actually carry a density
    // jump. Deep in either fluid rho is uniform, the buoyancy term is zero and
    // must stay zero; -snGrad(rho*gh) there is -rho*snGrad(gh), the
    // hydrostatic gradient that the p_rgh splitting has already removed. It is
    // a gradient, so the pressure absorbs it inside the region, but it drops
    // abruptly to zero at the edge of the region, and that step drives an
    // instability which took the standing wave to 10^37 m/s in 13 s.
    const scalarField& rhoI = rho.primitiveField();
    const scalar rhoRange = max(gMax(rhoI) - gMin(rhoI), SMALL);
    const surfaceScalarField snRho(fvc::snGrad(rho));
    const surfaceScalarField threshold
    (
        drhoTol_*rhoRange*mesh.deltaCoeffs()
    );

    scalarField& FI = F.primitiveFieldRef();
    const scalarField& snRhoI = snRho.primitiveField();
    const scalarField& thrI = threshold.primitiveField();
    const scalarField& snRhoGhI = snRhoGh.primitiveField();
    const scalarField& asgI = ghAssigned_.primitiveField();
    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();

    // reach = 1 only where the whole gradient stencil of the cell is
    // interface-based
    scalarField reach(asgI);

    forAll(own, facei)
    {
        const scalar asgOwn = asgI[own[facei]];
        const scalar asgNei = asgI[nei[facei]];
        reach[own[facei]] = min(reach[own[facei]], asgNei);
        reach[nei[facei]] = min(reach[nei[facei]], asgOwn);
    }

    scalarList nbrAsg;
    syncTools::swapBoundaryCellList(mesh, asgI, nbrAsg);

    forAll(mesh.boundary(), patchi)
    {
        const fvPatch& fvp = mesh.boundary()[patchi];
        const labelUList& fc = fvp.faceCells();
        const label offset = fvp.start() - nInt;

        if (fvp.coupled())
        {
            forAll(fc, i)
            {
                reach[fc[i]] = min(reach[fc[i]], nbrAsg[offset + i]);
            }
        }
        else
        {
            forAll(fc, i)
            {
                reach[fc[i]] = 0;
            }
        }
    }

    forAll(FI, facei)
    {
        if
        (
            reach[own[facei]] > 0.5
         && reach[nei[facei]] > 0.5
         && mag(snRhoI[facei]) > thrI[facei]
        )
        {
            FI[facei] = -snRhoGhI[facei];
        }
    }

    scalarList nbrReach;
    syncTools::swapBoundaryCellList(mesh, reach, nbrReach);

    surfaceScalarField::Boundary& Fbf = F.boundaryFieldRef();

    forAll(Fbf, patchi)
    {
        const fvPatch& fvp = mesh.boundary()[patchi];

        if (!fvp.coupled())
        {
            continue;
        }

        const fvsPatchScalarField& snRhoGhP = snRhoGh.boundaryField()[patchi];
        const fvsPatchScalarField& snRhoP = snRho.boundaryField()[patchi];
        const fvsPatchScalarField& thrP = threshold.boundaryField()[patchi];
        const labelUList& fc = fvp.faceCells();
        const label offset = fvp.start() - nInt;

        forAll(Fbf[patchi], i)
        {
            if
            (
                reach[fc[i]] > 0.5
             && nbrReach[offset + i] > 0.5
             && mag(snRhoP[i]) > thrP[i]
            )
            {
                Fbf[patchi][i] = -snRhoGhP[i];
            }
        }
    }
}


Foam::tmp<Foam::surfaceScalarField> Foam::densityWeightedGravityModel::force()
{
    calculateGhf();
    const volScalarField& rho = mesh_.lookupObject<volScalarField>("rho");

    tmp<surfaceScalarField> tF(-ghf_*fvc::snGrad(rho));

    if (densityWeighted_)
    {
        densityWeightedForce(rho, tF.ref());
    }

    return tF;
}


Foam::tmp<Foam::surfaceScalarField> Foam::densityWeightedGravityModel::ghf() const
{
    return ghf_;
}


Foam::tmp<Foam::volScalarField> Foam::densityWeightedGravityModel::gh() const
{
    return gh_;
}

// ************************************************************************* //
