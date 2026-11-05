/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
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
    fluxinator

Group
    grpMultiphaseSolvers

Description
    Two-phase Volume-of-Fluid solver with geometric interface advection
    (isoAdvector) and a consistent gravity force and momentum convection.

    The source is that of interIsoGoam, with these defaults instead of those
    of interIsoFoam:

    - gravityModel geometric, when there is no constant/gravityProperties
    - PIMPLE explicitConvection yes
    - PIMPLE ddtCorr no

    An entry in the case overrides a default. The convection scheme is set in
    system/fvSchemes, as in any solver. Explicit convection requires the
    interfaceUpwind scheme, and the solver warns without it. For waves use

        div(rhoPhi,U) Gauss laxWendroff U phi
                      interfaceUpwind alpha.water 1e-6 limitedLinearV 1;

\*---------------------------------------------------------------------------*/

#define FLUXINATOR

#include "../interIsoGoam/interIsoGoam.C"

// ************************************************************************* //
