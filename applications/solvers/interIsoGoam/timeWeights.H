// Time weights of the momentum convection and the gravity force (PIMPLE
// entries), evaluated once per time step:
//
//   convectionTheta: theta*div(rhoPhi, U^{n+1}) + (1 - theta)*div(rhoPhi, U^n)
//   gravityTheta:    theta*F^{n+1} + (1 - theta)*F^n (see pEqn.H)
//
// An entry in the PIMPLE dictionary is always used. The defaults depend on
// nOuterCorrectors and explicitConvection (default no, in fluxinator yes,
// see solverDefaults.H):
//
//   nOuter 1, explicitConvection no     convectionTheta 1,   gravityTheta 1
//   nOuter 1, explicitConvection yes    convectionTheta 0,   gravityTheta 1
//   nOuter > 1                          convectionTheta 1/2, gravityTheta 1/2
//
// The first line is interIsoFoam. The last line, with outerCorrBlending 1/2
// (alpha solver dictionary, default), is time-centred: alpha is advected with (phi^n + phi^{n+1})/2,
// and convection and gravity are averages of the start and end of the time
// step. With one outer corrector, U^{n+1} and F^{n+1} are not available
// before the pressure solution, and the convection term is explicit.
//
// The Lax-Wendroff face value (laxWendroff scheme) cancels the
// anti-diffusion (dt/2)(u & grad)^2 of forward Euler; the theta scheme has
// (1 - 2 theta) times this anti-diffusion, so the scheme scales its
// correction by the registered factor laxWendroffFactor = 1 - 2 theta
// (1 for explicit convection, 0 for the centred scheme).

const bool explicitConvection =
    pimple.dict().getOrDefault<bool>
    (
        "explicitConvection",
        defaultExplicitConvection
    );

const bool centredOuter = pimple.nCorrPIMPLE() > 1;

const scalar convTheta =
    pimple.dict().getOrDefault<scalar>
    (
        "convectionTheta",
        centredOuter ? 0.5 : (explicitConvection ? 0 : 1)
    );

const scalar gravityTheta =
    pimple.dict().getOrDefault<scalar>
    (
        "gravityTheta",
        centredOuter ? 0.5 : 1
    );

{
    if (!mesh.foundObject<uniformDimensionedScalarField>("laxWendroffFactor"))
    {
        auto* ptr = new uniformDimensionedScalarField
        (
            IOobject
            (
                "laxWendroffFactor",
                runTime.constant(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::REGISTER
            ),
            dimensionedScalar(dimless, 1)
        );
        ptr->store();
    }

    mesh.lookupObjectRef<uniformDimensionedScalarField>("laxWendroffFactor")
        .value() = 1 - 2*convTheta;

    static scalar convThetaPrinted = -1, gravityThetaPrinted = -1;
    if (convTheta != convThetaPrinted || gravityTheta != gravityThetaPrinted)
    {
        Info<< "Time weights: convectionTheta " << convTheta
            << ", gravityTheta " << gravityTheta
            << ", laxWendroffFactor " << 1 - 2*convTheta << nl << endl;
        convThetaPrinted = convTheta;
        gravityThetaPrinted = gravityTheta;
    }

    // Explicit convection amplifies velocity differences at the interface
    // unless the convected face velocity there is taken from the upwind cell
    static bool convectionSchemeChecked = false;
    if (convTheta == 0 && !convectionSchemeChecked)
    {
        convectionSchemeChecked = true;
        bool hasInterfaceUpwind = false;
        for (const token& tok : mesh.divScheme("div(rhoPhi,U)"))
        {
            hasInterfaceUpwind =
                hasInterfaceUpwind
             || (tok.isWord() && tok.wordToken() == "interfaceUpwind");
        }
        if (!hasInterfaceUpwind)
        {
            WarningInFunction
                << "Explicit convection without the interfaceUpwind scheme in"
                << " div(rhoPhi,U): velocity differences are amplified in"
                << " interface cells. Use, for example," << nl
                << "    div(rhoPhi,U) Gauss laxWendroff U phi interfaceUpwind"
                << " alpha.water 1e-6 limitedLinearV 1;" << nl << endl;
        }
    }
}
