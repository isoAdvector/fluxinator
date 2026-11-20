/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
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
    interIsoGoam

Group
    grpMultiphaseSolvers

Description
    The solver of Parts I and II, kept to reproduce their results. It is
    built from the source of fluxinator (../fluxinator), with the defaults
    of interIsoFoam instead of those of fluxinator:

    - interfaceForces.gravity.model basic (system/fvSchemes)
    - PIMPLE convectionNewTimeWeight 1 (implicit convection)
    - PIMPLE ddtCorr yes
    - nInitialReconstructions 0 (solver dictionary of alpha)

    With these defaults and one outer corrector it is interIsoFoam. An entry
    in the case overrides a default, so with the same settings interIsoGoam
    and fluxinator give identical results. The settings are listed in
    ../fluxinator/fluxinator.C.

    References:
    \verbatim
        Roenby, J., Missios, K., Møller, K., Aliyar, S. and Scheufler, H. (2026).
        Consistent interfacial flow numerics - Part I: The gravity force term.
        arXiv:2611.24841, https://arxiv.org/abs/2611.24841

        Roenby, J. (2026).
        Consistent interfacial flow numerics - Part II: The momentum
        convection term.
        arXiv:2611.24866, https://arxiv.org/abs/2611.24866
    \endverbatim

\*---------------------------------------------------------------------------*/

#define INTERISOGOAM

#include "../fluxinator/fluxinator.C"

// ************************************************************************* //
