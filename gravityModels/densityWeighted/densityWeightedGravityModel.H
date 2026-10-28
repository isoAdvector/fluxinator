/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2025 Henning Scheufler
    Copyright (C) 2026 Johan Roenby
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

Class
    Foam::densityWeightedGravityModel

Description
    Experimental variant of the geometric gravity model, kept for the record.

    Where the geometric model applies the buoyancy term as gh*snGrad(rho),
    this one replaces it by -snGrad(rho*gh) on the faces that can carry it.
    That makes the force a discrete gradient, so its loop imbalance vanishes
    and the spurious velocities of a moving, sloping interface fall by five to
    nine times, but it is not the same term: it differs by rho*grad(gh), the
    hydrostatic gradient that the p_rgh splitting has already removed. The
    price is a cell-scale sawtooth on the interface, which no choice of the
    options below removes, because a force that is both exactly gh*grad(rho)
    and exactly a gradient would need grad(gh) parallel to grad(rho), and at an
    interface they are perpendicular.

    See run/standingWave/experiments/FINDINGS.md. Not recommended for
    production use; the geometric model is the one the paper reports.

    Gravity force calculated based on reconstructed interface centres
    interpolated to mesh faces.

    g*h is set from the reconstructed interface in the interface cells and in
    nLayers rings of point neighbours around them; elsewhere it is g*h at the
    cell centre. The first ring uses a weighted average of the neighbouring
    interface centres. Each further ring uses the average of the values
    already assigned to its neighbours.

    With uncorrected snGrad schemes one layer is sufficient. With corrected
    (non-orthogonal) schemes, fvc::snGrad(rho) also reaches the faces between
    the first and second ring, so two layers are needed for an exact
    hydrostatic balance of a flat interface on non-orthogonal meshes.

Usage
    In constant/gravityProperties:
    \verbatim
    gravityModel    densityWeighted;
    nLayers         2;              // Optional, default 1
    ghAtCellCentre  false;          // Optional, default false
    densityWeighted false;          // Optional, default false
    drhoTol         1e-3;           // Optional, default 1e-3
    \endverbatim

    \c ghAtCellCentre samples the interface height above the cell centre rather
    than at the interface segment centroid, so that the linear interpolation to
    a face estimates the height at the face rather than mixing heights sampled
    at different horizontal positions. \c densityWeighted replaces the force
    -ghf*snGrad(rho) by -snGrad(rho*gh) on the faces whose two cells both carry
    an interface-based g*h. Both are
    aimed at interfaces that are not level; with a level interface g*h is
    constant near the interface and neither changes anything.

SourceFiles
    densityWeightedGravityModel.C

\*---------------------------------------------------------------------------*/

#ifndef densityWeightedGravityModel_H
#define densityWeightedGravityModel_H

#include "typeInfo.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "autoPtr.H"
#include "uniformDimensionedFields.H"
#include "gravityModel.H"

#include "zoneDistribute.H"
#include "markInterfaceRegion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class densityWeightedGravityModel Declaration
\*---------------------------------------------------------------------------*/

class densityWeightedGravityModel
:
    public gravityModel
{

        // g*h at faces based on geometrically reconstructed interface height
        surfaceScalarField ghf_;

        // g*h at cell centres based on geometrically reconstructed interface height
        volScalarField gh_;

        //- 1 where gh_ holds an interface-based value, 0 where it is the
        //  cell-centre height
        volScalarField ghAssigned_;

        //- Height reference field
        uniformDimensionedScalarField hRef_;

        //- Number of rings of cells around the interface cells that get the
        //  interface-based g*h
        label nLayers_;

        //- Take g*h of an interface cell as the interface height above the
        //  cell centre (the reconstructed plane projected along gravity)
        //  instead of the height of the interface segment centroid
        bool ghAtCellCentre_;

        //- Take the force to be -snGrad(rho*gh) rather than
        //  -ghf*snGrad(rho) on the faces that can carry it
        bool densityWeighted_;

        //- Density jump below which the force is left untouched, relative to
        //  the density range of the mesh
        scalar drhoTol_;

    // Private Member Functions

        //- Disallow copy construct
        densityWeightedGravityModel(const densityWeightedGravityModel&);

        //- Disallow default bitwise assignment
        void operator=(const densityWeightedGravityModel&);

        protected:

        //- Calculate the interface height
        void calculateGhf();

        //- Replace the force by -snGrad(rho*gh) near the interface
        void densityWeightedForce
        (
            const volScalarField& rho,
            surfaceScalarField& F
        ) const;


public:

    //- Runtime type information
    TypeName("densityWeighted");


    // Constructors

        //- Construct from components
        densityWeightedGravityModel
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~densityWeightedGravityModel()
    {}

    //- Calculate and return force on faces
    virtual tmp<surfaceScalarField> force();

    //- Return ghf
    tmp<surfaceScalarField> ghf() const;

    //- Return gh
    tmp<volScalarField> gh() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
