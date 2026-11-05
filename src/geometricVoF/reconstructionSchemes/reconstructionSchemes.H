/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2020 DLR
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
    Foam::reconstructionSchemes

Description
    Original code supplied by Henning Scheufler, DLR (2019)

SourceFiles
    reconstructionSchemes.C

\*---------------------------------------------------------------------------*/

#ifndef reconstructionSchemes_H
#define reconstructionSchemes_H

#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "dimensionedScalar.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "IOdictionary.H"
#include "DynamicField.H"
#include "MeshedSurface.H"
#include "MeshedSurfacesFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class reconstructionSchemes Declaration
\*---------------------------------------------------------------------------*/

class reconstructionSchemes
:
    public IOdictionary
{
    // Private Data

        //- Coefficients dictionary
        dictionary reconstructionSchemesCoeffs_;


protected:

    // Protected Data

        //- Reference to the VoF Field
        volScalarField& alpha1_;

        //- Reference to the face fluxes
        const surfaceScalarField& phi_;

        //- Reference to the velocity field
        const volVectorField& U_;

        //- Interface area normals
        volVectorField normal_;

        //- Interface centres
        volVectorField centre_;

        //- Is interface cell?
        boolList interfaceCell_;

        //- List of cell labels that have an interface
        DynamicField<label> interfaceLabels_;

        //- Store timeindex/iteration to avoid multiple reconstruction.
        //  Initialised to (-1, 0) so that alreadyReconstructed() returns
        //  false on the very first call, enabling lazy initialisation.
        mutable Pair<label> timeIndexAndIter_;

        //- Copies of interfaceCell and interfaceLabels captured by
        //  storeOld(), used together with normal_.oldTime() and
        //  centre_.oldTime() to cheaply reset the interface data via
        //  resetToOld() on later outer (PIMPLE) iterations within the
        //  same time step, without a full reconstruction. boolList and
        //  DynamicField have no built-in old-time mechanism of their
        //  own, unlike normal_ and centre_, so these must be stored
        //  explicitly.
        boolList interfaceCellOld_;
        DynamicField<label> interfaceLabelsOld_;


    // Protected Member Functions

        //- Is the interface already up-to-date?
        bool alreadyReconstructed(bool forceUpdate = true) const;

        //- Trigger reconstruct() if surface data has never been computed
        void ensureValid() const
        {
            if (timeIndexAndIter_.first() == -1)
                const_cast<reconstructionSchemes*>(this)->reconstruct(false);
        }

        //- No copy construct
        reconstructionSchemes(const reconstructionSchemes&) = delete;

        //- No copy assignment
        void operator=(const reconstructionSchemes&) = delete;


public:

    // Public Classes

        class interface
        :
            public meshedSurface
        {
            //- For every face the original cell in mesh
            labelList meshCells_;

            public:

                interface
                (
                    const pointField& pointLst,
                    const faceList& faceLst,
                    const labelList& meshCells
                )
                :
                    meshedSurface(pointLst, faceLst, surfZoneList()),
                    meshCells_(meshCells)
                {}

                interface
                (
                    pointField&& pointLst,
                    faceList&& faceLst,
                    labelList&& meshCells
                )
                :
                    meshedSurface
                    (
                        std::move(pointLst),
                        std::move(faceLst),
                        surfZoneList()
                    ),
                    meshCells_(std::move(meshCells))
                {}

                //- For every face, the original cell in mesh
                const labelList& meshCells() const
                {
                    return meshCells_;
                }
        };


public:

    //- Runtime type information
    TypeName("reconstructionSchemes");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            reconstructionSchemes,
            components,
            (
                volScalarField& alpha1,
                const surfaceScalarField& phi,
                const volVectorField& U,
                const dictionary& dict
            ),
            (alpha1, phi, U, dict)
        );


    // Selectors

        //- Return a reference to the selected phaseChange model
        static autoPtr<reconstructionSchemes> New
        (
            volScalarField& alpha1,
            const surfaceScalarField& phi,
            const volVectorField& U,
            const dictionary& dict
        );


    // Constructors

        //- Construct from components
        reconstructionSchemes
        (
            const word& type,
            volScalarField& alpha1,
            const surfaceScalarField& phi,
            const volVectorField& U,
            const dictionary& dict
        );


    //- Destructor
    virtual ~reconstructionSchemes() = default;


    // Member Functions

        //- Access to the model dictionary
        dictionary& modelDict() noexcept
        {
            return reconstructionSchemesCoeffs_;
        }

        //- Const access to the model dictionary
        const dictionary& modelDict() const noexcept
        {
            return reconstructionSchemesCoeffs_;
        }

        //- Reconstruct the interface
        virtual void reconstruct(bool forceUpdate = true) = 0;

        //- Record the current interface data (normal, centre,
        //  interfaceCell, interfaceLabels) as the old-time values for
        //  cheap restoration via resetToOld() on later outer (PIMPLE)
        //  iterations within this time step, avoiding the cost of a
        //  full reconstruction.
        void storeOld();

        //- Reset the interface data to the values captured by the most
        //  recent call to storeOld(), instead of reconstructing.
        void resetToOld();

        //- const-Reference to interface normals (triggers lazy init)
        const volVectorField& normal() const
        {
            ensureValid();
            return normal_;
        }

        //- const-Reference to interface centres (triggers lazy init)
        const volVectorField& centre() const
        {
            ensureValid();
            return centre_;
        }

        //- Reference to interface normals (for writing; does not trigger init)
        volVectorField& normal() noexcept
        {
            return normal_;
        }

        //- Reference to interface centres (for writing; does not trigger init)
        volVectorField& centre() noexcept
        {
            return centre_;
        }

        //- Does the cell contain interface (triggers lazy init)
        const boolList& interfaceCell() const
        {
            ensureValid();
            return interfaceCell_;
        }

        //- List of cells with an interface (triggers lazy init)
        const DynamicField<label>& interfaceLabels() const
        {
            ensureValid();
            return interfaceLabels_;
        }

        //- Map VoF Field in case of refinement
        virtual void mapAlphaField() const
        {}

        //- Generated interface surface points/faces.
        //  \note the points are disconnected!
        interface surface();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
