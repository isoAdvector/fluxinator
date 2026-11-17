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

Class
    Foam::geometricGravityModel

Description
    Gravity force calculated based on reconstructed interface centres
    interpolated to mesh faces.

    g*h is set from the reconstructed interface in the interface cells and in
    nLayers rings of point neighbours around them; elsewhere it is g*h at the
    cell centre. The first layer uses a weighted average of the neighbouring
    interface centres. Each further layer uses the average of the values
    already assigned to its neighbours.

    An interface lying exactly on mesh faces (volume fractions of exactly 0
    and 1, no interface cells) is handled by treating every face between two
    non-interface cells with densities on either side of the mean as an
    interface segment: both cells get the face centre and normal and count as
    interface cells before the layers are formed.

    With uncorrected snGrad schemes one layer is sufficient. With corrected
    (non-orthogonal) schemes, fvc::snGrad(rho) also reaches the faces between
    the first and second layer, so two layers are needed for an exact
    hydrostatic balance of a flat interface on non-orthogonal meshes.

Usage
    In system/fvSchemes:
    \verbatim
    interfaceForces
    {
        gravity
        {
            model           geometric;
            nLayers         2;  // Optional, default 2
        }
    }
    \endverbatim

SourceFiles
    geometricGravityModel.C

\*---------------------------------------------------------------------------*/

#ifndef geometricGravityModel_H
#define geometricGravityModel_H

#include "gravityModel.H"

#include "zoneDistribute.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class geometricGravityModel Declaration
\*---------------------------------------------------------------------------*/

class geometricGravityModel
:
    public gravityModel
{
    // Private Data

        //- g*h at faces based on the reconstructed interface height
        surfaceScalarField ghf_;

        //- g*h at cell centres based on the reconstructed interface height
        volScalarField gh_;

        //- Number of layers of cells around the interface cells that get the
        //  interface-based g*h
        const label nLayers_;


    // Private Member Functions

        //- Add the faces that are themselves the interface (volume fractions
        //  exactly 0 and 1 on either side) as interface segments of both
        //  their cells
        void addInterfaceFaces
        (
            boolList& interfaceCell,
            volVectorField& interfaceCentre,
            volVectorField& interfaceNormal
        ) const;

        //- Set gh_ in the interface cells and in layer 1, and mark them in
        //  assigned
        void setInterfaceCellsAndFirstLayer
        (
            const volVectorField& interfaceCentre,
            const volVectorField& interfaceNormal,
            const labelList& cellLayer,
            const boolList& nearInterface,
            zoneDistribute& exchangeFields,
            const scalar ghRef,
            volScalarField& assigned
        );

        //- Set gh_ in layers 2 to nLayers_: the average of the values already
        //  assigned to the point neighbours, one layer at a time
        void setFurtherLayers
        (
            const labelList& cellLayer,
            const boolList& nearInterface,
            zoneDistribute& exchangeFields,
            volScalarField& assigned
        );

        //- Calculate gh_ and ghf_ from the reconstructed interface
        void calculateGhf();


public:

    //- Runtime type information
    TypeName("geometric");


    // Constructors

        //- Construct from components
        geometricGravityModel
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~geometricGravityModel() = default;


    //- Calculate and return force on faces
    virtual tmp<surfaceScalarField> force();

    //- Return ghf
    tmp<surfaceScalarField> ghf() const;

    //- Return gh
    tmp<volScalarField> gh() const;

    //- Number of layers around the interface cells with an interface-based gh
    label nLayers() const
    {
        return nLayers_;
    }


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
