# Where the remaining spurious currents come from

All runs: interIsoGoam, geometric model, corrected schemes unless stated, one
time step of 0.01 s unless stated. The wave starts at rest at maximum
elevation, so the exact max|U| after one step is **0.0154 m/s**.

## 1. The case understates the model: it is solver-tolerance limited

With a **flat** interface on this mesh the exact solution is rest, so every
velocity is spurious. With `nLayers 2` and corrected schemes:

| p_rgh solver | max\|U\| after 20 steps |
|---|---|
| `tolerance 1e-07; relTol 0.05` (as shipped) | 6.9e-5 |
| `tolerance 1e-10; relTol 0` | 6.3e-8 |
| `tolerance 1e-12; relTol 0` | 6.0e-9 |

The 1e-5 floor was the linear solver, not the gravity model. With a tight
solver the model holds this refined mesh at rest to 6e-9 m/s, as on the meshes
of Section 5.

## 2. Two rings are needed here, as in Section 5

Flat interface, tight solver, 20 steps:

| setting | max\|U\| |
|---|---|
| basic, corrected | 2.96 |
| basic, uncorrected | 1.65 |
| geometric nLayers 1, corrected | 1.15 |
| geometric nLayers 2 or 3, corrected | 6.0e-9 |
| geometric nLayers 1 or 2, uncorrected | 2.4e-8 |

`snGradRhoReach` on this mesh: uncorrected reaches rings 0 and 1 only;
corrected reaches **4 cells in ring 2**, none beyond. Exactly the mechanism of
Section 5.4. `isoAlpha` and `plicRDF` give identical results to all digits.

## 3. In the wave, the water is right and only the air is wrong

Median over 12 interface offsets spanning one coarse cell, one time step:

| mesh | air max\|U\| median | air max | pure-water max\|U\| | water error |
|---|---|---|---|---|
| uniform 50x50 | 0.119 | 0.199 | 0.0152 | 1.3 % |
| uniform 100x100 | 0.141 | 0.200 | 0.0154 | 0.0 % |
| uniform 200x200 | 0.173 | 0.201 | 0.0159 | 3.2 % |
| refined band, interface inside it | 0.131 | 0.200 | 0.0154 | 0.3 % |
| refined square, interface crosses it | **0.680** | **1.499** | 0.0155 | 0.5 % |

The water solution is accurate to a few per cent everywhere, including on the
crossing mesh. All of the spurious velocity is in the air.

## 4. The baseline air error comes from the interface SLOPE, not its curvature

It is already there on a uniform mesh and does **not** converge: the worst air
velocity is 0.199, 0.200, 0.201 m/s for h = 0.02, 0.01, 0.005 m, so the residual
force per unit volume is independent of the cell size, as in Section 2.3.

A **tilted plane** initial condition, which has zero curvature everywhere,
produces the same kind of error, growing with the slope (first-step air max|U|,
median over 6 interface offsets):

| interface slope | uniform 50x50 | refined, crossing |
|---|---|---|
| 0 (flat) | 0.0000 | 0.0000 |
| 0.02 | 0.0088 | 0.125 |
| 0.05 | 0.0241 | 0.312 |
| 0.111 | 0.0787 | 0.325 |
| 0.2 | 0.1606 | 1.250 |

At slope 0.111, the slope the standing wave has where it crosses the
transition, the zero-curvature plane gives 0.079 / 0.325 against 0.119 / 0.680
for the wave. Curvature therefore accounts for at most a factor of two; the
error is first order in the slope and vanishes for a flat interface.

Edge case: with the plane exactly on a mesh face (z = 0.5) there are **no**
interface cells at all, the geometric model has nothing to reconstruct from, it
falls back to cell-centre heights everywhere, and the air velocity jumps to
4.5 m/s. Rare, but there is no graceful degradation.

## 5. The refinement transition only hurts where the interface crosses it

- refined band (transitions far from the interface): 1.10x the uniform error
- refined region below the interface: 1.00x
- refined square (interface crosses the transition): **5.7x median, 7.5x worst**

and the penalty scales with the interface slope at the crossing. With a 1 m
wavelength, so that the transitions at x = 0.25 and 0.75 fall either on the wave
extrema or on the nodes:

| transition falls on | amplification over uniform |
|---|---|
| wave extrema (interface flat there) | 1.4x |
| wave nodes (steepest slope) | 2.7x |

## 6. No available setting removes it

Four interface offsets, first-step air max\|U\| on the crossing mesh:

| setting | median |
|---|---|
| nLayers 1 | 1.018 |
| nLayers 2 | 0.670 |
| nLayers 3 | 0.670 |
| nLayers 4 | 0.670 |
| nLayers 2, uncorrected | 0.624 |
| nLayers 2, nCorrectors 5, nNonOrthogonalCorrectors 5 | 0.670 |

nLayers 2, 3 and 4 are bit-identical, so the reach is not what limits this.
The first-step velocity is exactly proportional to the time step (halving dt
halves max\|U\| to five digits), so it is an unbalanced force, not accumulation.

## 7. Interpretation: linear interpolation of gh is the wrong face value

Let G_i be the local interface height the model stores in cell i. With the fluid
at rest under a sloping interface the exact modified pressure is
psi_i = -rho_i G_i, which is a perfectly good cell field. For the discrete force
-gh_f (rho_N - rho_P)/|d_f| to be the discrete gradient of psi, the face value
must be

    gh_f = (rho_N G_N - rho_P G_P) / (rho_N - rho_P),

whereas `geometricGravityModel.C:332` uses `ghf_ = fvc::interpolate(gh_)`, the
linear interpolation w G_P + (1-w) G_N. Writing G_N = G_P + dG, the difference
between the two leaves an unbalanced force per unit volume of

    rho_f * dG / |d_f|,

that is, the interpolated density times the face gradient of the interface
height. It is zero for a flat interface (dG = 0), first order in the interface
slope, and independent of curvature - exactly what Section 4 measures. It also
explains the refinement transition: across a 2:1 face the two cells sample the
interface a larger distance apart, so dG is larger for the same slope, and it
explains why extra rings do nothing, since the problem is the *value* of gh and
not its *reach*.

The formula suggests a concrete remedy: form the face value with the
density-weighted expression above instead of interpolating, guarded where
rho_N - rho_P is small (there the density gradient vanishes and the value does
not matter). If it works it should make *any* interface whose height is
representable as a cell field exactly balanced, not just a flat one.

## 8. Practical note

The error varies by about a factor of ten with where the interface falls
relative to the cell faces, the same lattice sensitivity as in Section 3.3.
Single runs are not a reliable comparison; use medians over an ensemble of
interface offsets.

---

# Comparison with interFoam and interIsoFoam

All three solvers on the same meshes, tight pressure tolerance, median over 12
interface offsets. The domain maximum of |U| is in the air in every run.

## First time step: the gravity term alone

The interface has not moved yet, so this isolates the gravity discretisation.
The exact value is 0.0154 m/s.

| mesh | interFoam | interIsoFoam | interIsoGoam | improvement |
|---|---|---|---|---|
| uniform 50x50 | 0.773 | 0.773 | **0.119** | 6.5x |
| refined, interface crosses | 3.719 | 3.719 | **0.680** | 5.5x |

interFoam and interIsoFoam agree to four digits, as in the tilted box: they
share the gravity term, and the advection scheme cannot matter yet. The
refinement transition costs every solver a factor of about five.

## Twenty time steps (0.2 s), refined mesh

| solver | air peak | air at 0.2 s | water at 0.2 s | water error |
|---|---|---|---|---|
| interFoam | 3.719 | 0.729 | 0.2496 | 1.5 % |
| interIsoFoam | 3.914 | 1.920 | 0.2542 | 0.3 % |
| interIsoGoam | **1.425** | **0.505** | 0.2535 | **0.1 %** |

## One wave period (1.2 s), refined mesh, median over 6 offsets

| solver | peak air | mean air | interface band [cells] |
|---|---|---|---|
| interFoam | 3.562 | 0.841 | 412 |
| interIsoFoam | 4.180 | 2.266 | 104 |
| interIsoGoam | **1.872** | **0.763** | **88** |

interIsoGoam vs interIsoFoam: 2.2x lower peak, 3.0x lower mean.
interIsoGoam vs interFoam: 1.9x lower peak, but only 1.1x lower mean.

## How to read this

Against interIsoFoam, the same advection scheme with the standard gravity term,
the improvement is unambiguous: 5.5x at the first step and 2 to 3x sustained.

Against interFoam the first-step improvement is just as large, and the peak is
halved, but the *time-averaged* air velocity is almost the same. interFoam
reaches that by smearing: after one period its interface band is 412 cells
against 88 for interIsoGoam, and its volume fraction goes slightly unbounded
(min alpha = -2e-5). The smearing lowers the density contrast that drives the
spurious currents, which is the same effect reported for the tilted box in
Section 3, where interFoam's currents also grow more slowly than
interIsoFoam's. interIsoGoam achieves a lower peak *and* keeps the sharpest
interface and the most accurate water velocity.

---

# Density-weighted face value: implemented and tested

Two options were added to `gravityModels/geometric`, both defaulting to the
present behaviour, so nothing else changes unless they are switched on:

```
gravityModel    geometric;
nLayers         2;
ghAtCellCentre  false;   // sample g*h above the cell centre, not at the
                         // interface segment centroid
densityWeighted false;   // use (rho_N gh_N - rho_P gh_P)/(rho_N - rho_P)
drhoTol         1e-3;    // below this density jump, keep linear interpolation
```

Regression with the defaults: the flat interface still gives 6.0e-9, the
standing wave still gives 0.3209 at the first step, and the tilted box of
Section 3 still gives 4.02e-10. Unchanged to all digits.

## `ghAtCellCentre` does nothing

First-step air max|U| on the tilted plane, slope 0.111: 0.0787 with it off,
0.0788 with it on. For a plane cutting a square cell the interface segment
centroid already sits almost vertically above the cell centre, so the
projection moves the sample point hardly at all. The sampling location is not
what limits the accuracy.

## `densityWeighted` works, but only with the uncorrected gradient

Tilted plane, slope 0.111, refined mesh with the interface crossing the
transition, first-step air max|U|, median over 6 offsets:

| schemes | face value | air max\|U\| |
|---|---|---|
| corrected | linear | 3.25e-1 |
| corrected | density weighted | 7.75e-1 |
| corrected | density weighted, drhoTol 0.5 | 3.77e-1 |
| uncorrected | linear | 2.84e-1 |
| **uncorrected** | **density weighted** | **3.76e-2** |

The derivation assumes the two-point face gradient. With the corrected scheme
the density-weighted value is an extrapolation that can lie far outside
[gh_P, gh_N], and multiplying the non-orthogonal correction term by it makes
matters worse. With the uncorrected scheme it does what the derivation
predicts.

## On the standing wave it removes most of the error and all of the transition penalty

First-step air max\|U\| on the crossing mesh, median (worst) over 6 offsets,
physical value 0.0154 m/s:

| configuration | air max\|U\| |
|---|---|
| corrected, linear (as shipped) | 0.679 (1.499) |
| uncorrected, linear | 0.627 (1.346) |
| **uncorrected, density weighted** | **0.028 (0.030)** |

24x better on the median, 51x on the worst case, and 1.8x the physical value
rather than 44x. Over 20 time steps:

| configuration | air peak | air at 0.2 s | water at 0.2 s |
|---|---|---|---|
| crossing mesh, corrected + linear | 1.642 | 0.687 | 0.253 |
| crossing mesh, uncorrected + densityWeighted | **0.294** | **0.294** | 0.265 |
| uniform mesh, uncorrected + densityWeighted | 0.292 | 0.292 | 0.259 |

The air velocity on the crossing mesh (0.294) is now the same as on a uniform
mesh (0.292): **the refinement transition penalty is gone**, and the air
velocity is within 16 % of the physical 0.253 m/s instead of 171 % above it.
A flat interface still balances to 1.5e-9, so nothing is lost there.

## Caveat

This is only tested with uncorrected schemes. Section 5 needs the corrected
schemes on genuinely non-orthogonal meshes, and there the density-weighted
value currently hurts. A form that is consistent with the non-orthogonal
correction, or a limiter keeping gh_f within [gh_P, gh_N], would be the next
thing to try.

---

# Reconciling the density-weighted value with the corrected schemes

## The general form

Requiring that `-ghf*snGrad(rho)` equal `-snGrad(rho*gh)` with the *same* face
gradient operator fixes the face value to

    ghf = snGrad(rho*gh) / snGrad(rho)

With the uncorrected scheme this is the two-point expression of the previous
section, and with a level interface it reduces to the constant, so the linear
interpolation is recovered. `densityWeighted` now uses this quotient, guarded
where the density jump is small (`drhoTol`) and where either cell lacks an
interface-based g*h.

It does what it was meant to. Tilted plane, slope 0.111, crossing mesh,
first-step air max|U|:

| schemes | face value | air max\|U\| |
|---|---|---|
| corrected | linear | 3.25e-1 |
| corrected | two-point density weighted | 7.75e-1 |
| **corrected** | **snGrad quotient** | **3.75e-2** |
| uncorrected | snGrad quotient | 3.76e-2 |

The corrected schemes now match the uncorrected ones, an 8.6x gain over the
baseline, and the standing wave with the *corrected* schemes improves from
0.679 to 0.060 at the first step (11x) and from an air peak of 1.64 to 0.33
over 20 steps (5x), with the flat interface still balancing to 3.5e-9.

## But it must stay off for the Section 5 cases

| case | linear | snGrad quotient |
|---|---|---|
| tiltedBox, nLayers 1, 1 step | 4.02e-10 | 1.83e-10 |
| tiltedBoxRefined, nLayers 2, hydrostatic | 6.45e-09 | 6.35e-08 |
| skewBox, nLayers 1, hydrostatic | 7.68e-10 | 1.05e-01 |
| skewBox, nLayers 2 or 3, hydrostatic | 7.68e-10 | 2.63e-02 |

The skewed box regresses by seven orders of magnitude. The residual sits in the
cells next to the walls at x = 0 and x = 1, where the interface meets the
boundary. With a level interface the quotient equals the constant only if
`rho*gh` is proportional to `rho` over the whole stencil the correction uses,
and the Gauss cell gradient in a wall-adjacent cell also draws on the boundary
value of `gh`, which is the height of the face centre, not the interface
height. Extra rings move the internal part of the problem away but cannot
change the boundary contribution, which is why the error saturates at 2.6e-2.
Forcing the boundary value of `gh` to the cell value makes it far worse (4.5
m/s), so that is not the remedy either.

## Where this leaves it

`densityWeighted` is a large gain for a moving, sloping interface on meshes
that are orthogonal or nearly so, which is the standing wave case, and it must
be left off for the strongly non-orthogonal static cases of Section 5. It is
therefore an option, not a new default. Closing the gap needs a treatment of
the interface height on boundary faces, and of the cell gradient stencil of the
correction, that keeps `rho*gh` proportional to `rho` wherever the correction
reaches - a bigger change than the face value alone.


# Revision of densityWeighted (after the 20 s runs in run/standingWave/densWei)

**Outcome: abandoned.** The interface roughness below is not acceptable, and
the term is conceptually wrong - buoyancy is gh*grad(rho), not grad(gh*rho).
The `geometric` model has been restored to its plain form (nLayers only) and
this work is archived as a separate `densityWeighted` gravity model, which is
not recommended for use. What follows is the record of why.

Three observations from those runs: much lower spurious air velocities, a
visibly wigglier interface, and no spurious velocities at the boundaries. (The
`corrupted size vs. prev_size` at the end of every log, including the `basic`
one, is a duplicate-symbol clash between `libnewGeometricVoF` and OpenFOAM's
own `libgeometricVoF`, unrelated to this model and fixed elsewhere.)

## The interface roughness is real

Interface height per column over the uniformly refined band, at t = 20 s, as
the rms of its second difference along x:

| | rms(d2h) | max abs(d2h) |
|---|---|---|
| densityWeighted off | 1.93e-04 | 4.88e-04 |
| snGrad quotient | 1.88e-03 | 5.18e-03 |

Ten times rougher, and the profile is a cell-scale sawtooth along the whole
interface, not something localised at the refinement transition.

## The quotient was hiding a much simpler force

The quotient snGrad(rho*gh)/snGrad(rho) is multiplied by snGrad(rho) again, so
what is actually applied is -snGrad(rho*gh). On 99% of the faces where it was
used the quotient fell outside the interval spanned by the two cell values of
g*h, by up to 880 times the span, because rho_N - rho_P vanishes between two
interface cells of similar alpha while (rho*gh)_N - (rho*gh)_P does not.
Forming the force directly reproduces the quotient to 14 digits and removes
the singularity, the threshold it needed, and the unbounded ghf.

The force is now

    F = -snGrad(rho*gh)   on the faces that qualify
    F = -ghf*snGrad(rho)  elsewhere

and, being the face gradient of a cell field, it has zero loop imbalance by
construction.

## Two restrictions, each for a specific failure

**The gradient stencil must stay inside the interface region.** A corrected
snGrad adds a non-orthogonal correction built from the Gauss cell gradient,
which samples rho*gh on every face of the cell, and next to a physical boundary
on the boundary value of g*h - the face-centre height, not the interface
height. That, and nothing else, was the seven-order-of-magnitude skewed-box
regression: with an *uncorrected* snGrad the substitution is exact there
(3.27e-10 against 2.85e-10). So a cell only qualifies if it and all of its face
neighbours carry an interface-based g*h, and a cell against a physical boundary
never does.

**The face must carry a density jump.** Deep in either fluid rho is uniform,
the buoyancy term is zero and must stay zero, while -snGrad(rho*gh) there is
-rho*snGrad(gh), the hydrostatic gradient the p_rgh splitting has already
removed. It is a gradient, so the pressure absorbs it inside the region, but it
steps to zero at the region's edge. Without the `drhoTol` guard the standing
wave grew steadily from the first step and reached 1e37 m/s at t = 13.1 s.

## Where it now stands

Section 5, hydrostatic, corrected schemes:

| case | dw off | dw on |
|---|---|---|
| skewBox, nLayers 1 | 2.81e-10 | 2.91e-10 |
| skewBox, nLayers 2 | 2.81e-10 | 3.12e-10 |
| tiltedBoxRefined, nLayers 2 | 4.49e-09 | 6.75e-09 |
| tiltedBox, nLayers 1 | 4.02e-10 | 2.05e-10 |

The regression is gone: `densityWeighted` no longer has to be off for the
Section 5 cases.

Standing wave, corrected schemes, nLayers 2:

| | first step | peak over 20 steps | peak over 20 s | mean over last 5 s | rms(d2h) |
|---|---|---|---|---|---|
| dw off | 0.321 | 1.651 | 1.918 | 0.854 | 1.93e-04 |
| old quotient | 0.037 | 0.297 | 0.842 | 0.365 | 1.88e-03 |
| revised | 0.037 | 0.297 | 0.989 | 0.376 | 1.45e-03 |

nLayers 1 is now too thin for the stencil restriction and is unstable (peak
6.35); nLayers 3 and 4 give 1.73e-03 and 1.80e-03. nLayers 2 is the best.

## The roughness is intrinsic, not a bug

It only falls from 1.88e-03 to 1.45e-03 and does not improve with a wider band,
so it is not an edge effect. The substituted force differs from the buoyancy
term by -rho*grad(G), where G is the interface-height field. Along a sloping
interface the prefactor is the local density of the interface cells, which
swings between the two fluid densities from one cell to the next as alpha
varies, so the extra force is smooth in G but noisy in rho - which is the
sawtooth. It cannot be removed within this formulation: a force that is both
exactly -G*grad(rho) and exactly a gradient requires grad(G) parallel to
grad(rho), and at an interface they are perpendicular. The loop imbalance
proportional to the interface slope is irreducible; `densityWeighted` trades it
for a gradient force that is wrong by -rho*grad(G).

So it remains a trade: five to nine times less spurious velocity against a
1.4% cell-scale ripple on the interface (the ripple is 3.6e-04 m against a wave
amplitude of 0.025 m). It is still an option rather than a default, but it is
no longer one that breaks anything.
