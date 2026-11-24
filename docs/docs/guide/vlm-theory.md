---
sidebar_position: 1
---

# VLM Theory

The vortex lattice method solves inviscid, incompressible, irrotational flow
about thin lifting surfaces. Each surface is replaced by a sheet of horseshoe
vortices whose strengths are chosen so that the flow is tangent to the surface
at a set of control points. The result is a linear system in the circulations,
and everything downstream — loads, derivatives, trim, adjoints — follows from
solving it once.

This sheet states what Zephyr actually assembles and solves. The modelling
consequences are collected in [Validity and Limits](./limitations.md).

## Discretization

A lifting surface becomes a structured grid of trapezoidal panels: `nx`
chordwise by `ny` spanwise. Each spanwise column of panels is a **strip**, and
the strip is the unit the far field and the load reports are organized around.

Within each panel:

- the **bound vortex** lies on the panel quarter chord,
- two **trailing legs** run semi-infinitely downstream from the bound
  endpoints,
- the **control point** sits at the panel three-quarter chord,
- the **force point** is the bound-vortex centre.

Chordwise spacing follows AVL's CSPACER construction, so the quarter/three-
quarter placement is exact for each spacing law rather than approximated after
the fact. The section lift-slope factor (`CLAF`) moves the control point
chordwise, which is how it changes the effective section lift slope.

Camber enters as a slope, not as thickness: the panel normal is rotated by the
local camber slope of the airfoil at that chordwise station, plus the section
twist and any surface incidence. Thickness is carried on the geometry for
visualization and never affects the solve.

## The horseshoe velocity

The velocity a unit-strength horseshoe induces at a field point is the
Biot–Savart integral over its three legs. For the bound segment with endpoints
$\mathbf{r}_1$ and $\mathbf{r}_2$ measured from the field point,

$$
\mathbf{v}_{\text{bound}} = \frac{\Gamma}{4\pi}
\frac{(\mathbf{r}_1 \times \mathbf{r}_2)\,(r_1 + r_2)}
     {r_1 r_2\,(r_1 r_2 + \mathbf{r}_1 \cdot \mathbf{r}_2)}
$$

and the two trailing legs contribute their own closed-form terms. Zephyr
evaluates all three in `solver::horseshoe_vortex_velocity`, which is also the
kernel the gradient module differentiates.

Three details of that kernel are worth knowing, because they are what make the
derivatives and the symmetry behave:

**Core regularization.** A field point on a filament gives an unbounded
velocity. Each evaluation carries a core radius $\epsilon$ that enters as
$\sqrt{d^4 + \epsilon^4}$ in the denominators, which is smooth and keeps the
kernel differentiable everywhere. The radius is derived from panel size, not
from an absolute length, so it scales with the model. A target on the same
component as the source sees $10^{-4}$ of the source width; a target on another
component sees the larger of $10^{-4}$ of the source chord and $10^{-4}$ of the
source width.

**Scale invariance.** The degeneracy tests are dimensionless ratios —
$\sin^2\theta > \varepsilon$ for collinearity rather than a comparison of a
length-to-the-fourth quantity against a fixed epsilon. A model expressed in
millimetres gets the same lattice as the same model in metres.

**Reflection equivariance.** The two trailing legs are summed with each other
before being added to the bound term, and signed zeros are produced by
subtraction rather than unary negation. Floating-point addition is commutative
but not associative, so this is what makes a mirrored panel evaluate
bit-identically to its reflection — which is in turn what lets the mirror-reuse
and symmetric-factorization accelerations be exact rather than approximate.

## Compressibility

Subsonic compressibility is a Prandtl–Glauert correction applied inside the
kernel. With $\beta = \sqrt{1 - M_\infty^2}$, streamwise separations are scaled
by $1/\beta$ before the Biot–Savart evaluation, and the streamwise velocity
component is divided by $\beta$ again. Setting `FlightCondition::mach` to zero
recovers the incompressible kernel exactly.

Mach is a solver-level quantity: changing it invalidates the influence matrices
and their factorization. Mach derivatives are available analytically, and are
assembled from the same kernel rather than by re-solving at a perturbed Mach.

## The linear system

Enforcing flow tangency at every control point gives

$$
\sum_{j=1}^{N} A_{ij}\,\Gamma_j = -\,\mathbf{V}_i \cdot \hat{\mathbf{n}}_i
$$

where $A_{ij}$ is the normalwash the $j$-th unit horseshoe induces along the
$i$-th control point's normal, and $\mathbf{V}_i$ is the onset velocity there:
the unit freestream plus the rotational contribution of the body rates about the
reference point.

Zephyr keeps **two** influence operators, and the distinction is load-bearing:

| Operator | Points | Bound legs | Used for |
| --- | --- | --- | --- |
| Normalwash, $N \times N$ | control points | all, including the panel's own | the circulation solve |
| Velocity, $3N \times N$ | bound-force points | the panel's own leg excluded | force velocities |

Including a panel's own bound leg at its own control point is correct for the
boundary condition and wrong for the Kutta–Joukowski force, where the
self-induced contribution is singular and physically absent. Older versions of
this code used one operator for both. Both operators are assembled in the same
pass, and the adjoint differentiates both.

A surface flagged `NOALBE` (`set_affected_by_alpha_beta(false)`) gets a zero
right-hand side; it still induces velocity on everything else and still carries
loads, but the freestream does not drive it.

## Symmetry

There are three distinct mechanisms, and they are not interchangeable:

- **`y_duplicate`** on a surface generates a real mirrored copy with its own
  panels and unknowns. The lattice is the full aircraft.
- **`SymmetryConfig`** adds image planes during influence assembly. The images
  contribute velocity but have no unknowns of their own.
- **`SolverMode::HalfWing`** in `analysis::Analysis` keeps only one semispan's
  unknowns and folds the reflected influence into the system matrix, solving
  $A_{\text{same}} + A_{\text{reflected}}$. See
  [Selective Analysis](./selective-analysis.md).

Geometry and strip-load adjoints support explicit duplication but not image
planes: mirror a wing with `y_duplicate` if you intend to differentiate it.

## Loads

Two independent integrations are performed, and both are reported:

**Near field.** The Kutta–Joukowski theorem at each bound-force point,

$$
\mathbf{F}_i = \rho\,\Gamma_i\,(\mathbf{V}_i \times \boldsymbol{\ell}_i)
$$

with $\mathbf{V}_i$ the local velocity from the velocity operator and
$\boldsymbol{\ell}_i$ the bound-vortex vector. This gives all six coefficients
and the full strip-level distribution.

**Far field.** A Trefftz-plane integration over the trailing wake, which gives
lift, side force, and induced drag from the spanwise circulation distribution
alone:

$$
C_{D_i} = \frac{1}{S_{\text{ref}}}\sum_s \Gamma_s\,(w_{y,s}\,\Delta z_s - w_{z,s}\,\Delta y_s)
$$

The far-field induced drag is the more trustworthy of the two on a coarse
lattice, because it does not depend on resolving the chordwise pressure
distribution. [Loads and Coefficients](./loads.md) covers when the two
disagree and which to believe.

## Why this structure matters for gradients

The circulation solve is the only implicit step. Everything else — normals from
geometry, the right-hand side from the flight condition, forces from
circulation and velocity — is explicit. So a derivative of any scalar output
decomposes into an explicit term at fixed circulation plus one transposed solve
of the same matrix the primal used:

$$
\frac{dJ}{dx} = \frac{\partial J}{\partial x}\bigg|_\Gamma
+ \boldsymbol{\lambda}^{\!\top}\!\left(\frac{\partial \mathbf{b}}{\partial x} - \frac{\partial A}{\partial x}\boldsymbol{\Gamma}\right),
\qquad
A^{\!\top}\boldsymbol{\lambda} = \frac{\partial J}{\partial \boldsymbol{\Gamma}}
$$

The factorization is already in hand from the primal solve, so one adjoint costs
a triangular solve rather than a factorization — and a hundred geometric design
variables cost one adjoint, not a hundred re-solves. That is the subject of
[Analytical and Adjoint Gradients](./automatic-differentiation.md).
