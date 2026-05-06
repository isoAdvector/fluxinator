/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2025 Henning Scheufler
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

#include "gravityReconModel.H"
#include "gravityMeshObject.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"
#include "reconstructionSchemes.H"
#include "cutFaceAdvect.H"
#include "processorPolyPatch.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(gravityReconModel, 0);
    addToRunTimeSelectionTable(gravityModel, gravityReconModel, components);
}


Foam::vector Foam::gravityReconModel::closestDist(const point p, const vector n ,const vector centre)
{
    vector normal = n/mag(n);
    return p - normal*((p - centre) & normal);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::gravityReconModel::gravityReconModel
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
    )
{
    calculateGhf();
}


// * * * * * * * * * * * * * * Public Access Member Functions  * * * * * * * //

void Foam::gravityReconModel::calculateGhf()
{

    const fvMesh& mesh = gh_.mesh();
    zoneDistribute& exchangeFields = zoneDistribute::New(mesh);

    // Read only if mesh changed would be clever
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
        reconstructionSchemes& surf = mesh.lookupObjectRef<reconstructionSchemes>("reconstructionScheme");

        surf.reconstruct(false);

        const volVectorField& faceCentre = surf.centre();
        const volVectorField& faceNormal = surf.normal();

        boolList nextToInterface(mesh.nCells(),false);
        markInterfaceRegion markIF(mesh);

        markIF.markCellsNearSurf(surf.interfaceCell(),1,nextToInterface);

        exchangeFields.setUpCommforZone(nextToInterface,true);

        Map<vector > mapCentres =
            exchangeFields.getDatafromOtherProc(nextToInterface,faceCentre);
        Map<vector > mapNormal =
            exchangeFields.getDatafromOtherProc(nextToInterface,faceNormal);

        const labelListList& stencil = exchangeFields.getStencil();

        forAll(surf.interfaceCell(),celli)
        {
            if(mag(faceNormal[celli]) != 0)
            {
                vector closeP = closestDist(mesh.C()[celli],-faceNormal[celli],faceCentre[celli]);
                gh_[celli] = closeP & g.value();
            }
            else if(nextToInterface[celli])
            {
                // the to the interface
                vector averageSurfP = vector::zero;
                scalar avgWeight = 0;
                const point p = mesh.C()[celli];

                forAll(stencil[celli],i)
                {
                    const label& gblIdx = stencil[celli][i];
                    vector n = -exchangeFields.getValue(faceNormal,mapNormal,gblIdx);
                    if (mag(n) != 0)
                    {
                        n /= mag(n);
                        vector c =
                            exchangeFields.getValue(faceCentre,mapCentres,gblIdx);
                        vector distanceToIntSeg = (c - p);
                        vector closeP = closestDist(p,n,c);
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
                        averageSurfP += closeP * weight;
                        avgWeight += weight;
                    }
                }

                if (avgWeight != 0)
                {
                    averageSurfP /= avgWeight;
                    gh_[celli] = averageSurfP & g.value();
                }

            }
            else
            {
                // do nothing
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


Foam::tmp<Foam::surfaceScalarField> Foam::gravityReconModel::force()
{
    calculateGhf();
    const volScalarField& rho = mesh_.lookupObject<volScalarField>("rho");
    return -ghf_*fvc::snGrad(rho);
}

// ************************************************************************* //
