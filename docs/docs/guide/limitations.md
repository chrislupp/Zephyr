---
sidebar_position: 14
---

# Validity and Limits

A vortex-lattice solution is cheap because it assumes a great deal. This sheet
states what Zephyr assumes, what it has been checked against, and which pieces
of its own API are present but inert. Read it before trusting a number near the
edge of the model.

## What the method assumes

**Potential flow.** Inviscid, incompressible, irrotational. Viscosity enters
only as the additive profile-drag terms described in
[Loads and Coefficients](./loads.md), which are corrections applied after the
fact and are not coupled back into the circulation.

**Thin surfaces.** Lifting surfaces are camber sheets. Thickness is carried on
the geometry for visualization and never affects the solve, so thickness drag,
thickness-induced interference, and transonic area ruling are all outside the
model.

**Attached flow.** There is no stall, no separation, and no maximum lift. The
lift curve is linear in angle of attack as far as you care to extend it. A
converged solution at 25 degrees is still a converged solution; it is simply not
an aircraft.

**A fixed, planar wake.** The trailing legs run straight downstream from the
bound endpoints. There is no wake relaxation, no roll-up, and no force-free wake
iteration. Post-stall wakes, close formation flight, and strongly interacting
wakes are not modelled.

**Small-displacement boundary conditions.** Flow tangency is enforced on the
undeflected lattice with rotated normals. Control deflection rotates the
boundary condition, not the geometry, which is accurate for moderate
deflections and degrades as the deflection grows.

**Subsonic compressibility only.** A Prandtl–Glauert correction, and nothing
more. Mach at or above 1 is rejected outright rather than silently returning the
incompressible answer. There is no transonic, shock, or supersonic model.

**Steady, quasi-steady rates.** Body rates enter as a steady rotational onset
velocity. There is no unsteady aerodynamics: no reduced-frequency dependence, no
apparent-mass terms, no indicial response.

## Specific consequences

| If you want | The honest answer |
| --- | --- |
| $C_{L,\max}$ or stall angle | Not available. Nothing in the model saturates. |
| Profile or form drag | Only through `CDp` and `CDCL` polars you supply. |
| Fuselage aerodynamics | Not modelled. `Body` geometry is inert — see below. Represent what you need as a `NOALBE` lifting shell. |
| Hinge moments | Not computed. `ControlDerivatives::hinge_moment` is always zero. |
| Transonic drag rise | Outside the model. |
| Ground effect | Representable as a `SymmetryConfig` image plane in `z`, but note that image planes are unsupported by the geometry and strip-load adjoints. |
| Flutter or gust response | Not available; the aerodynamics are steady. |
| Near-field versus far-field drag | Prefer `far_field_induced_drag` on a coarse lattice. It needs no chordwise resolution. |

## Present but inert

These exist in the API and currently do nothing. They are listed here so you do
not assume otherwise:

| Item | Status |
| --- | --- |
| `Body`, `BodySection`, and the AVL `BODY`/`BSEC` blocks | parsed and stored; no contribution to the solution |
| `Surface::set_has_wake(false)` (`NOWAKE`) | recorded per panel; the solver does not act on it |
| `Surface::set_contributes_to_loads(false)` (`NOLOAD`) | recorded; the force integration does not exclude it |
| `ControlDerivatives::hinge_moment` | always zero |
| `StripLoads::induced_angle`, `StripLoads::induced_velocity` | never filled; the VTK `InducedAlpha` array is zero with them |
| `TrimPresets::steady_turn`'s bank argument | accepted for signature parity; no turn kinematic relation is imposed |
| AVL's `INDEX` keyword | not recognized; use `COMPONENT` |

Per-panel induced velocity *is* available, in `SolutionState::induced_velocity`.

## Restrictions inside the library

**Geometry and strip-load adjoints** reject two configurations rather than
returning an approximate derivative:

- lattices using `SymmetryConfig` image planes — duplicate explicitly instead;
- meshes carrying section drag polars, because the strip-area and quarter-chord
  tangents the viscous chain needs do not exist yet.

Flight-condition and control derivatives have neither restriction.

**Half-wing mode** requires surfaces duplicated about `y = 0`, zero sideslip,
zero roll and yaw rate, and symmetric control reflection signs; antisymmetric
derivatives are rejected. It never silently falls back to a full solve. The
complete list is in [Selective Analysis](./selective-analysis.md).

**Mesh-node and section gradients** hold reference quantities fixed. If your
reference area or chord depends on the design variables, that chain rule is
yours to supply.

**Trim** needs a nonsingular local Jacobian; bound-active control allocation is
not provided, and `MeshTrimSolver` reports a failure rather than approximating
one. There is no continuation warm start.

**Rolling and yawing moment coefficients** are normalized by the reference
chord, not the span. This differs from AVL — see
[Conventions](../getting-started/conventions.md).

## What has been checked

Aerodynamic results are compared against **independent OptVL 2.1.1 fixtures**,
committed in `tests/fixtures/avl_reference/` together with the AVL geometry
files they were produced from and the OptVL version and commit they came from.
The cases are:

| Case | Condition | Panels |
| --- | --- | --- |
| `rectangular` | $\alpha = 4^\circ$, $\beta = 2^\circ$, $M = 0.2$ | 8 |
| `rectangular_claf` | the same, with a section lift-slope factor | 8 |
| `rectangular_sine` | the same, with sine chordwise spacing | 8 |
| `wing_tail` | $\alpha = 5^\circ$, $M = 0.3$ | 44 |
| `ucrm_5_8` | $\alpha = 2^\circ$, symmetric half model | 40 |
| `ucrm_17_40` | $\alpha = 2^\circ$, symmetric half model | 680 |

Each case checks lift, drag, and pitching moment, their alpha derivatives, and
the full circulation vector, in both full and half-wing modes where both apply.
Coefficients are required to agree to better than `2e-8` absolute and
derivatives to `2e-7`.

Within the library, selective and legacy paths are compared against each other
with a tolerance of `1e-10 + 1e-8 * |reference|`, and analytical gradients are
compared against directional finite differences with separate allowances for
truncation and roundoff. Those bounds apply to the tested cases, not to
arbitrarily ill-conditioned systems; the selective algorithms may differ in
rounding from the full solver.

The frozen v0.2.0 oracle in `tests/fixtures/parity/` remains in use for span
geometry and dynamics, which the port did not affect. **Its aerodynamic values
contain port errors** and are not a reference for loads. The behavioral
differences that matter:

- Normalwash includes every bound vortex at the control points, while force
  velocities exclude the self-bound leg and are evaluated at the bound-force
  points. Moments use those same points. v0.2.0 used one operator for both.
- Section cosine spacing and `CLAF` follow AVL's CSPACER construction.
- Direct node meshes keep quarter/three-quarter panel placement and physically
  twisted nodes, so importing sections and supplying twisted nodes need not
  produce identical lattices.

A manifest in `tests/fixtures/parity/ctest-manifest.json` maps the non-Philote
v0.2.0 test cases to their Rust replacements.

No timing or adoption claims are published. If performance matters to your case,
measure it — see [Performance](./performance.md).

## Reporting a disagreement

If Zephyr disagrees with AVL or OptVL on a case you can share, the fixtures
above are the pattern to follow: the geometry file, the operating condition, and
the reference numbers. Check the normalization conventions first — reference
area, and the chord-versus-span moment normalization — because that is where
most apparent disagreements turn out to live.
