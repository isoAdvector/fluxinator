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
    convectionTest

Description
    Kinematic test of the momentum convection of a VoF solver with a moving
    density jump (Part II). A 2D velocity field with one empty direction n is
    prescribed,
        u(r, t) = U0 + a t + (Omega - beta |d|^2) n x d + S ((r - r0).e_s) e_f,
    with d = (r - c) in the plane and e_s = n x e_f: uniform translation and
    acceleration, solid-body rotation, a vortex with curvature in its velocity
    profile (beta = Omega2), shear, or any superposition. The face fluxes come
    from the stream function psi (u = grad(psi) x n), so they are divergence
    free to round-off on any mesh. alpha is advected with isoAdvector, and
    rhoPhi is isoAdvector's mass flux.

    Each variant in system/convectionTestDict solves the momentum equation
    without pressure,
        ddt(rho, Uv) + div(rhoPhi, Uv) = rho* a_exact,
    with a_exact = dU/dt + grad(u).u. Named as the time weights of fluxinator,
    forceNewTimeWeight w_f gives the force density rho* = w_f rho^{n+1}
    + (1 - w_f) rho^n (default 1), and convectionNewTimeWeight w_c the
    convection w_c div(rhoPhi, Uv^{n+1}) + (1 - w_c) div(rhoPhi, Uv^n), with
    any ddt and div scheme. Mode oneStep starts every step from the exact
    field (the error of one step); mode carried carries Uv from step to step
    (the accumulated error).

    Output, for e = Uv - u(reference, t^{n+1}): per step and variant
    postProcessing/convectionTest/<variant>.dat holds
        t maxAll maxWater maxInterface rmsInterface maxAir momentumError
    (water: alpha >= 0.99 at both time levels, air: <= 0.01 at both, interface:
    the rest; momentumError = |sum V (rho Uv - rho u)|). Error fields
    err_<variant> are written at write times.

    Reference (convectionTestDict entry reference, default cellCentre):
      cellCentre    u and a_exact at the cell centre. In interface cells the
                    exact mass-averaged velocity differs from it by
                    O(|grad u| dx) for rotation and shear.
      massCentroid  u and a_exact at the mass centroid of the cell, from a PLIC
                    cut with the normal grad(alpha) (Youngs), not from the
                    advector's reconstruction, so that the advection does not
                    depend on the reference. For a linear u this is the exact
                    mass-averaged cell velocity.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "isoAdvector.H"
#include "surfaceIteratorPLIC.H"
#include "ddtScheme.H"
#include "convectionScheme.H"
#include "OFstream.H"
#include "StringStream.H"
#include "emptyFvPatch.H"
#include "zeroField.H"
#include "uniformDimensionedFields.H"

using namespace Foam::fluxinatorVoF;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

struct prescribedFlow
{
    vector n, U0, a, c, r0, ef, es;
    scalar Omega, beta, S;

    prescribedFlow(const dictionary& d, const vector& nEmpty)
    :
        n(nEmpty),
        U0(d.getOrDefault<vector>("U0", Zero)),
        a(d.getOrDefault<vector>("acceleration", Zero)),
        c(d.getOrDefault<vector>("rotationCentre", Zero)),
        r0(d.getOrDefault<vector>("shearOrigin", Zero)),
        ef(d.getOrDefault<vector>("shearDirection", vector(1, 0, 0))),
        es(Zero),
        Omega(d.getOrDefault<scalar>("Omega", 0)),
        beta(d.getOrDefault<scalar>("Omega2", 0)),
        S(d.getOrDefault<scalar>("shearRate", 0))
    {
        ef -= (ef & n)*n;
        ef /= mag(ef);
        es = n ^ ef;
    }

    vector V(scalar t) const
    {
        return U0 + a*t;
    }

    vector u(const vector& r, scalar t) const
    {
        const vector rp = r - (r & n)*n;
        const vector cp = c - (c & n)*n;
        return
            V(t) - (V(t) & n)*n
          + (Omega - beta*magSqr(rp - cp))*(n ^ (rp - cp))
          + S*((r - r0) & es)*ef;
    }

    vector accel(const vector& r, scalar t) const
    {
        // dU/dt + grad(u).u, grad(u) by central differences
        const scalar h = 1e-6;
        const vector uc = u(r, t);
        vector conv(Zero);
        for (direction j = 0; j < 3; ++j)
        {
            vector e(Zero);
            e[j] = h;
            conv += uc[j]*(u(r + e, t) - u(r - e, t))/(2*h);
        }
        return a - (a & n)*n + conv;
    }

    scalar psi(const vector& r, scalar t) const
    {
        const vector rp = r - (r & n)*n;
        const vector cp = c - (c & n)*n;
        const scalar s = (r - r0) & es;
        return
            ((n ^ V(t)) & r)
          - 0.5*Omega*magSqr(rp - cp)
          + 0.25*beta*sqr(magSqr(rp - cp))
          + 0.5*S*sqr(s);
    }
};


void setFlux
(
    const prescribedFlow& flow,
    scalar t,
    surfaceScalarField& phi
)
{
    const fvMesh& mesh = phi.mesh();
    const pointField& pts = mesh.points();
    const faceList& faces = mesh.faces();

    auto faceFlux = [&](const label facei)
    {
        const face& f = faces[facei];
        scalar s = 0;
        forAll(f, i)
        {
            const point& p1 = pts[f[i]];
            const point& p2 = pts[f.nextLabel(i)];
            s += flow.psi(0.5*(p1 + p2), t)*(flow.n & (p2 - p1));
        }
        return s;
    };

    scalarField& phiIn = phi.primitiveFieldRef();
    forAll(phiIn, facei)
    {
        phiIn[facei] = faceFlux(facei);
    }
    auto& phibf = phi.boundaryFieldRef();
    forAll(phibf, patchi)
    {
        const label start = mesh.boundaryMesh()[patchi].start();
        forAll(phibf[patchi], i)
        {
            phibf[patchi][i] =
                isA<emptyFvPatch>(mesh.boundary()[patchi])
              ? 0
              : faceFlux(start + i);
        }
    }
}


void setVelocity
(
    const prescribedFlow& flow,
    scalar t,
    volVectorField& U
)
{
    const fvMesh& mesh = U.mesh();
    forAll(U, celli)
    {
        U[celli] = flow.u(mesh.C()[celli], t);
    }
    auto& Ubf = U.boundaryFieldRef();
    forAll(Ubf, patchi)
    {
        if
        (
            !Ubf[patchi].coupled()
         && !isA<emptyFvPatch>(mesh.boundary()[patchi])
        )
        {
            const vectorField& Cf = mesh.boundary()[patchi].Cf();
            vectorField v(Cf.size());
            forAll(v, i)
            {
                v[i] = flow.u(Cf[i], t);
            }
            Ubf[patchi] == v;
        }
    }
    U.correctBoundaryConditions();
}


//- Point per cell at which the exact velocity is evaluated: the cell centre,
//  or the mass centroid from a PLIC cut with the normal grad(alpha)
tmp<pointField> refPoints
(
    const bool massCentroid,
    const volScalarField& alpha1,
    const scalar rho1,
    const scalar rho2,
    surfaceIteratorPLIC& sIter
)
{
    const fvMesh& mesh = alpha1.mesh();
    auto tx = tmp<pointField>::New(mesh.C().primitiveField());
    if (!massCentroid)
    {
        return tx;
    }
    pointField& x = tx.ref();
    const volVectorField normal(fvc::grad(alpha1));
    const scalar tol = 1e-12;
    forAll(x, celli)
    {
        const scalar a = alpha1[celli];
        const vector& nw = normal[celli];
        if (a <= tol || a >= 1 - tol || mag(nw) < SMALL)
        {
            continue;
        }
        sIter.vofCutCell(celli, a, 1e-12, 100, nw);
        const scalar V = mesh.V()[celli];
        const scalar Vw = sIter.subCellVolume();
        const scalar Va = V - Vw;
        const point& xw = sIter.subCellCentre();
        const point xa = (V*mesh.C()[celli] - Vw*xw)/Va;
        x[celli] = (rho1*Vw*xw + rho2*Va*xa)/(rho1*Vw + rho2*Va);
    }
    return tx;
}

} // End namespace Foam


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Kinematic test of VoF momentum convection with a prescribed linear"
        " flow and isoAdvector (see the source header)"
    );

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    IOdictionary testDict
    (
        IOobject
        (
            "convectionTestDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const Vector<label>& gD = mesh.geometricD();
    vector nEmpty(Zero);
    for (direction d = 0; d < 3; ++d)
    {
        if (gD[d] == -1)
        {
            nEmpty[d] = 1;
        }
    }
    if (mag(nEmpty) != 1)
    {
        FatalErrorInFunction
            << "convectionTest needs a 2D mesh with one empty direction"
            << exit(FatalError);
    }

    const prescribedFlow flow(testDict.subDict("flow"), nEmpty);
    const word mode(testDict.getOrDefault<word>("mode", "oneStep"));
    if (mode != "oneStep" && mode != "carried")
    {
        FatalErrorInFunction
            << "mode " << mode << ", valid: oneStep carried"
            << exit(FatalError);
    }
    const word reference
    (
        testDict.getOrDefault<word>("reference", "cellCentre")
    );
    if (reference != "cellCentre" && reference != "massCentroid")
    {
        FatalErrorInFunction
            << "reference " << reference << ", valid: cellCentre massCentroid"
            << exit(FatalError);
    }
    const bool massCentroid = (reference == "massCentroid");

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
    const dimensionedScalar rho1
    (
        "rho1",
        dimDensity,
        transportProperties.subDict(phases[0]).get<scalar>("rho")
    );
    const dimensionedScalar rho2
    (
        "rho2",
        dimDensity,
        transportProperties.subDict(phases[1]).get<scalar>("rho")
    );

    volScalarField alpha1
    (
        IOobject
        (
            "alpha." + phases[0],
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    setVelocity(flow, runTime.value(), U);

    surfaceScalarField phi
    (
        IOobject
        (
            "phi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimVolume/dimTime, Zero)
    );
    phi.setOriented();
    setFlux(flow, runTime.value(), phi);
    Info<< "max |div(phi)| = " << gMax(mag(fvc::div(phi))().primitiveField())
        << " 1/s" << endl;

    volScalarField rho
    (
        IOobject("rho", runTime.timeName(), mesh),
        alpha1*rho1 + (1 - alpha1)*rho2
    );
    rho.oldTime().oldTime();

    // Variants
    const dictionary& variantsDict = testDict.subDict("variants");
    const wordList names(variantsDict.toc());
    PtrList<volVectorField> Uv(names.size());
    List<word> ddtSchemes(names.size()), divSchemes(names.size());
    scalarList convectionWeight(names.size()), forceWeight(names.size());
    PtrList<OFstream> files(names.size());
    mkDir(runTime.globalPath()/"postProcessing"/"convectionTest");

    forAll(names, vi)
    {
        const dictionary& vd = variantsDict.subDict(names[vi]);
        ddtSchemes[vi] = vd.getOrDefault<string>("ddt", "Euler");
        divSchemes[vi] = vd.get<string>("div");
        convectionWeight[vi] = vd.get<scalar>("convectionNewTimeWeight");
        forceWeight[vi] = vd.getOrDefault<scalar>("forceNewTimeWeight", 1);

        Uv.set
        (
            vi,
            new volVectorField
            (
                IOobject
                (
                    "U_" + names[vi],
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                U
            )
        );
        Uv[vi].oldTime().oldTime();

        files.set
        (
            vi,
            new OFstream
            (
                runTime.globalPath()/"postProcessing"/"convectionTest"
               /(names[vi] + ".dat")
            )
        );
        files[vi].precision(8);
        files[vi]
            << "# " << names[vi] << ": ddt " << ddtSchemes[vi]
            << ", div " << divSchemes[vi]
            << ", convectionNewTimeWeight " << convectionWeight[vi]
            << ", forceNewTimeWeight " << forceWeight[vi]
            << ", mode " << mode << ", reference " << reference << nl
            << "# t maxAll maxWater maxInterface rmsInterface maxAir"
            << " momentumError" << nl;
    }

    alpha1.oldTime();
    isoAdvector advector(alpha1, phi, U);

    // The convection time weight of the variant being solved, registered as in
    // fluxinator for the laxWendroff scheme
    uniformDimensionedScalarField convectionNewTimeWeight
    (
        IOobject
        (
            "convectionNewTimeWeight",
            runTime.constant(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::REGISTER
        ),
        dimensionedScalar(dimless, Zero)
    );
    zeroField Sp, Su;
    surfaceIteratorPLIC sIter(mesh, 1e-12);

    // Initial velocity of the variants at the reference points
    if (massCentroid)
    {
        const pointField x0
        (
            refPoints(true, alpha1, rho1.value(), rho2.value(), sIter)
        );
        forAll(Uv, vi)
        {
            forAll(x0, celli)
            {
                Uv[vi][celli] = flow.u(x0[celli], runTime.value());
            }
        }
    }

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        ++runTime;
        const scalar tNew = runTime.value();
        Info<< "Time = " << runTime.timeName() << endl;

        // Advect with the start-of-step flux and velocity
        advector.advect(Sp, Su);
        const surfaceScalarField rhoPhi(advector.getRhoPhi(rho1, rho2));
        rho == alpha1*rho1 + (1 - alpha1)*rho2;

        const volScalarField& rho0 = rho.oldTime();

        // Exact acceleration and exact new velocity at the reference points
        // (for the mass centroid: from alpha^{n+1})
        const pointField x1
        (
            refPoints(massCentroid, alpha1, rho1.value(), rho2.value(), sIter)
        );
        volVectorField aEx
        (
            IOobject("aEx", runTime.timeName(), mesh),
            mesh,
            dimensionedVector(dimVelocity/dimTime, Zero)
        );
        forAll(aEx, celli)
        {
            aEx[celli] = flow.accel(x1[celli], tNew);
        }
        volVectorField UexNew("UexNew", U);
        setVelocity(flow, tNew, UexNew);
        forAll(x1, celli)
        {
            UexNew[celli] = flow.u(x1[celli], tNew);
        }

        const scalarField& a1 = alpha1.primitiveField();
        const scalarField& a0 = alpha1.oldTime().primitiveField();

        forAll(names, vi)
        {
            volVectorField& W = Uv[vi];
            // Boundary values: exact at the new time
            auto& Wbf = W.boundaryFieldRef();
            forAll(Wbf, patchi)
            {
                if
                (
                    !Wbf[patchi].coupled()
                 && !isA<emptyFvPatch>(mesh.boundary()[patchi])
                )
                {
                    Wbf[patchi] == UexNew.boundaryField()[patchi];
                }
            }

            tmp<fv::ddtScheme<vector>> tddt
            (
                fv::ddtScheme<vector>::New
                (
                    mesh,
                    IStringStream(ddtSchemes[vi])()
                )
            );
            convectionNewTimeWeight.value() = convectionWeight[vi];
            tmp<fv::convectionScheme<vector>> tconv
            (
                fv::convectionScheme<vector>::New
                (
                    mesh,
                    rhoPhi,
                    IStringStream(divSchemes[vi])()
                )
            );

            const scalar wf = forceWeight[vi];
            const tmp<volScalarField> trhoStar
            (
                wf == 1 ? tmp<volScalarField>(rho)
              : wf == 0 ? tmp<volScalarField>(rho0)
              : wf*rho + (1 - wf)*rho0
            );
            const volVectorField F(trhoStar()*aEx);

            // As UEqn.H of fluxinator: the implicit operator only for w_c > 0
            const scalar wc = convectionWeight[vi];
            fvVectorMatrix WEqn
            (
                wc == 1
              ? tmp<fvVectorMatrix>
                (
                    tddt.ref().fvmDdt(rho, W) + tconv().fvmDiv(rhoPhi, W)
                 ==
                    F
                )
              : wc == 0
              ? tmp<fvVectorMatrix>
                (
                    tddt.ref().fvmDdt(rho, W)
                  + tconv().fvcDiv(rhoPhi, W.oldTime())
                 ==
                    F
                )
              : tmp<fvVectorMatrix>
                (
                    tddt.ref().fvmDdt(rho, W)
                  + wc*tconv().fvmDiv(rhoPhi, W)
                  + (1 - wc)*tconv().fvcDiv(rhoPhi, W.oldTime())
                 ==
                    F
                )
            );
            WEqn.solve();

            // Errors
            const vectorField e(W.primitiveField() - UexNew.primitiveField());
            scalar maxAll = 0, maxW = 0, maxI = 0, maxA = 0, sumI = 0;
            label nI = 0;
            vector mom(Zero);
            forAll(e, celli)
            {
                const scalar me = mag(e[celli]);
                maxAll = max(maxAll, me);
                const scalar amin = min(a1[celli], a0[celli]);
                const scalar amax = max(a1[celli], a0[celli]);
                if (amin >= 0.99)
                {
                    maxW = max(maxW, me);
                }
                else if (amax <= 0.01)
                {
                    maxA = max(maxA, me);
                }
                else
                {
                    maxI = max(maxI, me);
                    sumI += sqr(me);
                    ++nI;
                }
                mom += mesh.V()[celli]*rho[celli]*e[celli];
            }
            reduce(maxAll, maxOp<scalar>());
            reduce(maxW, maxOp<scalar>());
            reduce(maxI, maxOp<scalar>());
            reduce(maxA, maxOp<scalar>());
            reduce(sumI, sumOp<scalar>());
            reduce(nI, sumOp<label>());
            reduce(mom, sumOp<vector>());

            if (Pstream::master())
            {
                files[vi]
                    << tNew << ' ' << maxAll << ' ' << maxW << ' ' << maxI
                    << ' ' << (nI ? Foam::sqrt(sumI/nI) : 0) << ' ' << maxA
                    << ' ' << mag(mom) << nl;
            }

            if (runTime.writeTime())
            {
                volVectorField err("err_" + names[vi], W - UexNew);
                err.write();
            }

            if (mode == "oneStep")
            {
                W == UexNew;
            }
        }

        // Prescribed flow at the new time (start of the next step), at the
        // cell centres whatever the reference
        setVelocity(flow, tNew, U);
        setFlux(flow, tNew, phi);

        runTime.write();
        if (runTime.writeTime())
        {
            rho.write();
        }
        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
