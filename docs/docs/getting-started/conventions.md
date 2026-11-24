---
sidebar_position: 3
---

# Conventions

A vortex-lattice result is only as meaningful as the frame, the units, and the
signs it is reported in. This sheet states all three in one place. Everything
here is the library's actual behavior, not a recommendation.

## Axes

Geometry lives in a right-handed body frame:

| Axis | Direction |
| --- | --- |
| `x` | aft, along the fuselage |
| `y` | starboard |
| `z` | up |

Angular rates are the body-axis components of the same frame: roll about `x`,
pitch about `y`, yaw about `z`, each positive by the right-hand rule. Rates act
about the solver's reference point, which is the aircraft's moment reference
point unless you override it with `VortexLatticeSolver::set_reference_point`.

The freestream direction is a unit vector built from angle of attack and
sideslip:

$$
\hat{\mathbf{V}}_\infty =
\begin{pmatrix} \cos\alpha\cos\beta \\ \sin\beta \\ \sin\alpha\cos\beta \end{pmatrix}
$$

Because `x` is aft and `z` is up, positive $\alpha$ tilts the onset flow
upward across the wing, which is the sign that produces positive lift.

## Units

Zephyr is unit-agnostic and converts nothing. Pick one consistent system for
length, area, mass, density, and dynamic pressure and stay in it. Everything is
double precision.

One split is deliberate and is worth memorizing, because it is the one that
silently produces a wrong answer:

| Quantity | Unit |
| --- | --- |
| `Section` twist, `Surface` incidence offset | **degrees** |
| `alpha`, `beta` | **radians** |
| Body rates | dimensionless (below) |
| Control deflections | **radians** |
| Trim variables | **radians** |
| Every derivative with respect to an angle | **per radian** |
| Section twist derivatives | **per degree** |

Degrees survive on the geometry side because that is what AVL geometry files
carry, and reading a file should not change the numbers in it.
`core::DEGREES_TO_RADIANS` and `RADIANS_TO_DEGREES` are provided for the
conversions you will write anyway.

## Body rates

`FlightCondition::roll_rate`, `pitch_rate`, and `yaw_rate` are the components of
an angular-velocity vector that is added to a **unit** freestream when the
boundary condition is formed:

$$
\mathbf{V}_i = \hat{\mathbf{V}}_\infty + \boldsymbol{\omega} \times (\mathbf{r}_i - \mathbf{r}_{\text{ref}})
$$

So a rate of 1.0 means an angular velocity of one radian per the time the
freestream takes to cover one length unit of the geometry: $\Omega \ell / V_\infty$
with $\ell$ equal to one length unit. The rate derivatives in
`StabilityDerivatives` are per unit of that same quantity.

`EigenmodeAnalysis` converts them to a dimensional system with the conventional
half-span and half-chord factors, $b/2V$ laterally and $c/2V$ longitudinally, so
there it reads them as derivatives with respect to $\hat p = pb/2V$,
$\hat q = qc/2V$, and $\hat r = rb/2V$. The two readings agree only when the
relevant reference length is two geometry length units; check the scaling for
your own model before taking absolute mode frequencies from a hand-assembled
`StabilityDerivatives`.

## Reference quantities

Coefficients are normalized by what you supply, never by anything inferred from
the planform:

| Field | Used for |
| --- | --- |
| `reference_area` | every force and moment coefficient |
| `reference_chord` | **every** moment coefficient, including roll and yaw |
| `reference_span` | aspect ratio and span efficiency only |
| `moment_reference_point` | the point all moments are taken about |

The second row differs from AVL, which normalizes rolling and yawing moments by
span. Zephyr uses the reference chord for all three moment coefficients. When
comparing `Cl` or `Cn` against AVL output, scale by $c_{\text{ref}}/b_{\text{ref}}$.

`Aircraft` carries these for section-based geometry;
`geometry::ReferenceQuantities` carries the same four for a direct mesh, and
`ReferenceQuantities::from(&aircraft)` extracts them.

## Forces and moments

Near-field panel forces come from the Kutta–Joukowski theorem evaluated at the
bound-vortex force point, with the panel's own bound leg excluded from its
velocity. Strip moments are taken about the strip quarter chord and then
transferred to the moment reference point.

`AerodynamicResults::body_force` and `body_moment` are coefficient vectors in
body axes. `coefficients` is the conventional six:

| Field | Definition |
| --- | --- |
| `lift` | $C_L = C_{F_z}\cos\alpha - C_{F_x}\sin\alpha$ |
| `drag` | $C_D = C_{F_x}\cos\alpha + C_{F_z}\sin\alpha$ |
| `side` | $C_Y = C_{F_y}$ |
| `roll` | $C_l$, the **negated** right-hand moment about `x` |
| `pitch` | $C_m$, the right-hand moment about `y` |
| `yaw` | $C_n$, the **negated** right-hand moment about `z` |

Lift and drag are stability-axis quantities, obtained by rotating the body force
about `y` through $\alpha$ alone; sideslip is not rotated out. The two sign
flips on `roll` and `yaw` are the conventional aerodynamic signs: a positive
`roll` is right-wing-down, a positive `yaw` is nose-right.

Where the library reports **dimensional** loads — `trim::AeroQuantity::Fx`
through `Mz` and the Philote discipline's `Fx`–`Mz` outputs — it uses body axes
and mathematical right-hand moments, so `Mx` and `Mz` do *not* carry those sign
flips. External loads handed to the trim solver must use the same axes and the
same moment reference point.

## Panel and strip ordering

A lattice is a flat list of horseshoe vortices and a flat list of spanwise
strips. Within a surface, panels run chordwise inside each strip, and strips run
spanwise. Across the lattice, each input surface is followed immediately by its
mirrored copy when `y_duplicate` is set. `VortexLattice::surface_info` gives the
contiguous vortex and strip ranges for every generated surface, including
mirrors; `SurfaceMeshInfo::original_surface_index` points back to the surface a
mirror came from, indexing `surface_info` rather than `Aircraft::surfaces`.

That ordering is stable and is what every array in the library is indexed by:
circulation, induced velocity, strip loads, mesh-node gradients, and strip-load
cotangents.

## Compressibility

A subsonic Prandtl–Glauert correction is applied through the influence
coefficients, driven by `FlightCondition::mach`. There is no transonic or
supersonic model; see [Validity and Limits](../guide/limitations.md).
