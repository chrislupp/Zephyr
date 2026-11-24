---
sidebar_position: 5
---

# Loads and Coefficients

`ForceCalculator` turns a circulation into loads. It performs two independent
integrations — near field and far field — and reports both, because they fail in
different ways and the disagreement between them is diagnostic.

```rust
let results = ForceCalculator::new().compute_for_aircraft(
    aircraft, solver.mesh(), solution, condition,
);
```

`compute_for_aircraft` takes reference quantities and the constant profile-drag
coefficient from the `Aircraft`. `compute_with_references` is the direct-mesh
entry point and takes a `ReferenceQuantities` instead. Both return an
`AerodynamicResults`; an invalid solution or an empty lattice yields a default
result with `is_valid == false` rather than an error.

## Near field

Each panel force is the Kutta–Joukowski resultant at the bound-vortex force
point,

$$
\mathbf{F}_i = 2\,\Gamma_i\,(\mathbf{V}_i \times \boldsymbol{\ell}_i)
$$

in the solver's unit-freestream normalization, where $\mathbf{V}_i$ is the
freestream plus the induced velocity at that point and $\boldsymbol{\ell}_i$ is
the bound-vortex vector. The panel's own bound leg is excluded from
$\mathbf{V}_i$ — that is what the separate velocity operator exists for.

`ForceCalculator::set_include_induced_velocity(false)` drops the induced term
and evaluates the forces against the freestream alone. The default is to include
it.

Strip moments are taken about the strip quarter-chord point, then transferred to
the moment reference point by $\mathbf{r} \times \mathbf{F}$.

This path gives all six coefficients, per-strip and per-surface records, and
near-field induced drag. It is sensitive to chordwise resolution, because it
depends on where along the chord the load actually sits.

## Far field

The Trefftz-plane integration uses only the **total circulation of each strip**,
taken at the strip's last chordwise bound vortex, and the trailing-leg geometry.
It gives lift, side force, and induced drag:

$$
C_L = \frac{2}{S_{\text{ref}}}\sum_s \Gamma_s \Delta y_s,
\qquad
C_Y = -\frac{2}{S_{\text{ref}}}\sum_s \Gamma_s \Delta z_s,
\qquad
C_{D_i} = \frac{1}{S_{\text{ref}}}\sum_s \Gamma_s (w_{y,s}\Delta z_s - w_{z,s}\Delta y_s)
$$

and span efficiency

$$
e = \frac{C_L^2 + C_Y^2}{\pi\,\mathrm{AR}\,C_{D_i}},
\qquad \mathrm{AR} = \frac{b_{\text{ref}}^2}{S_{\text{ref}}}
$$

Because it needs no chordwise information, far-field induced drag converges with
far fewer chordwise panels than the near-field value does. **On a coarse
lattice, trust `far_field_induced_drag` over `near_field_induced_drag`.**

Strips on the same surface see no core radius in the Trefftz sum; strips on
different surfaces get one, taking the larger of the two scaled lengths. It
defaults to 0.01 of the source chord and 0.1 of the source width and is
configurable with `set_trefftz_core_radius(chord_fraction, width_fraction)`.
That is what keeps a wing and a closely trailing tail from producing a singular
far-field integral.

The regularization inside the influence kernel is separate and much tighter:
$10^{-4}$ of the source width within a component, and the larger of $10^{-4}$ of
the source chord and width across components.

## Reading the result

| Field | Meaning |
| --- | --- |
| `is_valid` | false means the inputs were rejected; every number is a default |
| `coefficients` | the conventional six, as defined in [Conventions](../getting-started/conventions.md) |
| `body_force`, `body_moment` | coefficient vectors in body axes |
| `near_field_induced_drag` | near-field `drag` less `viscous_drag` |
| `viscous_drag` | area-weighted strip polar contribution, plus the aircraft constant |
| `far_field_lift`, `far_field_side`, `far_field_induced_drag` | Trefftz values |
| `span_efficiency` | as above; zero when induced drag is negligible |
| `reference_*`, `moment_reference_point` | the normalization actually used |
| `surface_loads` | one record per **generated** surface, mirrors included |

`AerodynamicResults::aspect_ratio()` and `ideal_induced_drag()` are convenience
functions on the reference quantities; the latter is $C_L^2/\pi\mathrm{AR}$, the
elliptic ideal to compare against.

## Strip and surface records

`SurfaceLoads` holds a surface's coefficients, its dimensionless force and
moment vectors, its area, and a `StripLoads` for each of its strips.

`StripLoads` is where spanwise distributions come from:

| Field | Meaning |
| --- | --- |
| `chord`, `width`, `area` | strip geometry |
| `leading_edge_position`, `quarter_chord_position` | strip reference points |
| `coefficients` | the six, normalized by the **reference** area and chord |
| `body_force`, `body_moment` | per-strip-area coefficient vectors |
| `normal_span_loading` | load normal to the strip, the span-loading curve |
| `viscous_drag` | this strip's polar contribution |

`StripLoads` also declares `induced_angle` and `induced_velocity`. The force
integration does not currently fill them; they stay zero, and the VTK
`InducedAlpha` array is zero with them. Per-panel induced velocity is available
in `SolutionState::induced_velocity`, an `N × 3` matrix in body axes indexed by
panel.

Surfaces appear in lattice order, each original followed by its mirror, matching
`VortexLattice::surface_info`. Strips within a surface run spanwise.

A caution on normalization: strip coefficients divide by the **aircraft**
reference area and chord, not by the strip's own. Summing `strip.coefficients`
over a surface reproduces the surface coefficient; dividing by the strip's own
area to get a local $c_l$ is a conversion you have to do yourself.

## Viscous drag

Two additive mechanisms, both inviscid-result-plus-correction and neither
coupled back into the circulation:

- `Aircraft::set_profile_drag_coefficient` adds one constant to the aircraft
  drag coefficient.
- Section and surface `CDCL` polars give each strip a $C_d(C_l)$ evaluated at
  its own local lift coefficient. The strip contribution is applied along the
  freestream direction, so it rotates with $\alpha$, and it is clamped at zero
  from below.

The total appears as `viscous_drag` and is also folded into `coefficients.drag`.
`near_field_induced_drag` is the difference, so it is the quantity to compare
against `far_field_induced_drag`.

See [Geometry](./geometry.md) for the polar's quadratic fit and its deliberate
lack of clamping outside the fitted range.

## Dimensional loads

`AerodynamicResults` is entirely dimensionless. To get forces, multiply by
`dynamic_pressure * reference_area` — and by `reference_chord` as well for
moments, including roll and yaw.

Where the library produces dimensional loads itself — `trim::AeroQuantity::Fx`
through `Mz`, and the Philote discipline's `Fx`–`Mz` outputs — it does that
multiplication for you and reports mathematical right-hand moments in body
axes, without the conventional `Cl`/`Cn` sign flips. External loads supplied to
the trim solver must match that convention and reference point.

## Lower-level access

`ForceCalculator::trefftz_plane` runs the far-field integration alone and
returns `(lift, side, induced_drag, efficiency)`, which is useful when the
near-field pass is not wanted. For distributed loads differentiated against an
arbitrary functional — the aeroelastic case — see
[Coupled Interfaces](./coupled-interfaces.md).
