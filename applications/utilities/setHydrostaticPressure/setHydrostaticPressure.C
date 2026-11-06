/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

Application
    setHydrostaticPressure

Group
    grpPreProcessingUtilities

Description
    Sets p_rgh to the hydrostatic field of a flat interface at rest,

        p_rgh = -rho*((g & x_I) - ghRef),

    where rho = alpha1*rho1 + (1 - alpha1)*rho2 is formed from the current
    volume fraction field, x_I is a point on the interface, and ghRef is the
    reference geopotential formed from constant/hRef exactly as in the gravity
    models of interIsoGoam. The interface pressure is taken as zero.

    For a flat interface this field is an exact discrete hydrostatic state for
    any face-normal gradient scheme that is linear in the cell values, provided
    the same scheme is used for p_rgh and rho and the gravity model evaluates
    the geopotential at the interface on all faces where the density gradient
    is non-zero. Starting a simulation from it avoids the start-up currents
    that arise when p_rgh starts at zero, in particular with the non-orthogonal
    corrections of the corrected schemes, which are lagged by one solution.

    Gravity is read from constant/g unless -gravity is given. The interface
    position is taken from -height, from -origin, or from the entry origin in
    system/setAlphaFieldDict (the dictionary of setAlphaField, as used by the
    constantHeight gravity model). Only a flat interface perpendicular to
    gravity is a state of rest; the utility warns if it is given another one.

    The densities are read from the first two phases in
    constant/transportProperties unless -rho1 and -rho2 are given.

    The internal field is set from the volume fraction, and the boundary values
    follow it wherever the boundary condition does not fix the value, so that
    the written field is hydrostatic up to and including the boundaries. A
    boundary condition that fixes the value is part of the case and is kept.

Usage
    \b setHydrostaticPressure [OPTIONS]

    Examples:
    \verbatim
        setHydrostaticPressure
        setHydrostaticPressure -height 0.5
        setHydrostaticPressure -gravity '(0 0 -9.81)' -height 0.5
        setHydrostaticPressure -origin '(0 0 0.5)' -alpha alpha.water
    \endverbatim

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "timeSelector.H"
#include "uniformDimensionedFields.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Set p_rgh to the hydrostatic field of a flat interface at rest"
    );

    timeSelector::addOptions_singleTime();  // Single-time options

    argList::addOption
    (
        "gravity",
        "vector",
        "Gravitational acceleration, e.g. '(0 0 -9.81)'. Default: constant/g"
    );
    argList::addOption
    (
        "height",
        "scalar",
        "Interface elevation measured opposite to gravity"
    );
    argList::addOption
    (
        "origin",
        "point",
        "Point on the interface, e.g. '(0 0 0.5)'"
    );
    argList::addOption
    (
        "dict",
        "file",
        "Dictionary holding the entry origin. Default: system/setAlphaFieldDict"
    );
    argList::addOption
    (
        "alpha",
        "name",
        "Name of the volume fraction field. Default: alpha.water"
    );
    argList::addOption
    (
        "p",
        "name",
        "Name of the pressure field. Default: p_rgh"
    );
    argList::addOption("rho1", "scalar", "Density of the first phase");
    argList::addOption("rho2", "scalar", "Density of the second phase");

    #include "setRootCase.H"
    #include "createTime.H"

    // Set the time from -time / -latestTime / -constant if given, else the
    // start time of controlDict
    timeSelector::setTimeIfPresent(runTime, args);

    #include "createMesh.H"

    const word alphaName = args.getOrDefault<word>("alpha", "alpha.water");
    const word pName = args.getOrDefault<word>("p", "p_rgh");

    // Gravitational acceleration
    vector gValue(Zero);
    if (args.found("gravity"))
    {
        gValue = args.get<vector>("gravity");
    }
    else
    {
        uniformDimensionedVectorField g
        (
            IOobject
            (
                "g",
                runTime.constant(),
                runTime,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );
        gValue = g.value();
    }

    if (mag(gValue) < SMALL)
    {
        FatalErrorInFunction
            << "Gravity is zero: no hydrostatic pressure to set"
            << exit(FatalError);
    }

    // Reference geopotential, as in the gravity models
    uniformDimensionedScalarField hRef
    (
        IOobject
        (
            "hRef",
            runTime.constant(),
            runTime,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        dimensionedScalar(dimLength, Zero)
    );

    const scalar ghRef =
        (gValue & (cmptMag(gValue)/mag(gValue)))*hRef.value();

    // Point on the interface
    if (args.found("height") && args.found("origin"))
    {
        FatalErrorInFunction
            << "Specify either -height or -origin, not both"
            << exit(FatalError);
    }

    point interfacePoint(Zero);

    if (args.found("height"))
    {
        interfacePoint = -gValue/mag(gValue)*args.get<scalar>("height");
    }
    else if (args.found("origin"))
    {
        interfacePoint = args.get<point>("origin");
    }
    else
    {
        const fileName dictPath = args.getOrDefault<fileName>("dict", "");

        IOdictionary interfaceDict
        (
            dictPath.empty()
          ? IOobject
            (
                "setAlphaFieldDict",
                runTime.system(),
                runTime,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
          : IOobject
            (
                dictPath,
                runTime,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        interfacePoint = interfaceDict.get<point>("origin");

        vector normal(Zero);
        if (interfaceDict.readIfPresent("normal", normal) && mag(normal) > SMALL)
        {
            const scalar alignment =
                mag(normal & gValue)/(mag(normal)*mag(gValue));

            if (alignment < 1 - 1e-6)
            {
                WarningInFunction
                    << "The interface normal " << normal
                    << " is not parallel to gravity " << gValue << nl
                    << "    Such an interface is not a state of rest, and the"
                    << " field written here is not an equilibrium." << endl;
            }
        }
    }

    // Densities
    scalar rho1 = 0;
    scalar rho2 = 0;

    if (args.found("rho1") && args.found("rho2"))
    {
        rho1 = args.get<scalar>("rho1");
        rho2 = args.get<scalar>("rho2");
    }
    else
    {
        IOdictionary transportProperties
        (
            IOobject
            (
                "transportProperties",
                runTime.constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        const wordList phases(transportProperties.get<wordList>("phases"));

        if (phases.size() < 2)
        {
            FatalErrorInFunction
                << "Expected two phases in " << transportProperties.name()
                << ", found " << phases
                << exit(FatalError);
        }

        rho1 = dimensionedScalar
        (
            "rho",
            dimDensity,
            transportProperties.subDict(phases[0])
        ).value();

        rho2 = dimensionedScalar
        (
            "rho",
            dimDensity,
            transportProperties.subDict(phases[1])
        ).value();

        rho1 = args.getOrDefault<scalar>("rho1", rho1);
        rho2 = args.getOrDefault<scalar>("rho2", rho2);
    }

    const scalar ghInterface = (gValue & interfacePoint) - ghRef;

    Info<< "Gravity        : " << gValue << nl
        << "Interface point: " << interfacePoint << nl
        << "hRef           : " << hRef.value() << nl
        << "g & x_I - ghRef: " << ghInterface << nl
        << "Densities      : " << rho1 << ", " << rho2 << nl << endl;

    Info<< "Time = " << runTime.timeName() << nl << endl;

    volScalarField alpha1
    (
        IOobject
        (
            alphaName,
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    volScalarField p_rgh
    (
        IOobject
        (
            pName,
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    p_rgh.primitiveFieldRef() =
       -(
            alpha1.primitiveField()*rho1
          + (1 - alpha1.primitiveField())*rho2
        )*ghInterface;

    // Let the boundary values follow the new internal field, except where the
    // boundary condition fixes the value. Writing the field also writes the
    // boundary entries, and a stale value there (value uniform 0 on a
    // fixedFluxPressure patch, say) is not the hydrostatic state: in the
    // skewed box it drives a spurious velocity of 7e-3 m/s in ten time steps,
    // against 8e-10 m/s with boundary values taken from the adjacent cells.
    volScalarField::Boundary& p_rghBf = p_rgh.boundaryFieldRef();

    forAll(p_rghBf, patchi)
    {
        fvPatchScalarField& pf = p_rghBf[patchi];

        if (!pf.coupled() && !pf.fixesValue())
        {
            pf == pf.patchInternalField();
        }
    }

    Info<< "Writing " << p_rgh.name() << ": min = "
        << gMin(p_rgh.primitiveField()) << ", max = "
        << gMax(p_rgh.primitiveField()) << nl << endl;

    p_rgh.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
