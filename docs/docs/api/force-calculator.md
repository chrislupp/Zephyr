---
sidebar_position: 7
---

# ForceCalculator

`aero::ForceCalculator` turns a circulation into loads. It performs two
independent integrations — near field and far field — and reports both. See
[Loads and Coefficients](../guide/loads.md) for the formulae and for when the
two disagree.

```rust
let results = ForceCalculator::new()
    .compute_for_aircraft(aircraft, solver.mesh(), solution, condition);
```

## Entry points

| Method | Takes |
| --- | --- |
| `compute_for_aircraft` | an `Aircraft`; uses its reference quantities and its constant profile drag |
| `compute_with_references` | a `ReferenceQuantities`; the direct-mesh path, with no aircraft-level profile drag |
| `trefftz_plane` | the far-field integration alone, returning `(lift, side, induced_drag, efficiency)` |

An invalid solution or an empty lattice yields a default `AerodynamicResults`
with `is_valid == false` rather than an error. **Check `is_valid` before reading
the numbers**; a default result is all zeros, which is a plausible-looking lie.

## Settings

| Method | Effect |
| --- | --- |
| `set_include_induced_velocity(bool)` | include the induced term in the force velocity. Default `true`. |
| `set_trefftz_core_radius(chord_fraction, width_fraction)` | regularization between strips on **different** surfaces, taking the larger of the two scaled lengths. Defaults: 0.01 of the source chord and 0.1 of the source width. |

Strips on the same surface see no Trefftz core radius. The cross-surface core is
what keeps a wing and a closely trailing tail from producing a singular
far-field integral.

## `AerodynamicResults`

| Field | Meaning |
| --- | --- |
| `is_valid` | false means the inputs were rejected |
| `coefficients` | the conventional six |
| `body_force`, `body_moment` | coefficient vectors in body axes |
| `near_field_induced_drag` | near-field `drag` less `viscous_drag` |
| `viscous_drag` | area-weighted strip polar contribution plus the aircraft constant |
| `far_field_lift`, `far_field_side`, `far_field_induced_drag` | Trefftz values |
| `span_efficiency` | $(C_L^2 + C_Y^2)/(\pi\,\mathrm{AR}\,C_{D_i})$; zero when induced drag is negligible |
| `reference_*`, `moment_reference_point` | the normalization actually used |
| `surface_loads` | one record per **generated** surface, mirrors included |

`aspect_ratio()` is $b^2/S$ and `ideal_induced_drag()` is $C_L^2/\pi\mathrm{AR}$,
the elliptic ideal to compare against.

Prefer `far_field_induced_drag` on a coarse lattice: it needs no chordwise
resolution, while the near-field value depends on resolving where along the
chord the load sits.

## `SurfaceLoads` and `StripLoads`

Surfaces appear in lattice order, each original followed immediately by its
mirror. Strips run spanwise within a surface.

Strip coefficients divide by the **aircraft** reference area and chord, not by
the strip's own, so they sum to the surface coefficient. Getting a local $c_l$
is a conversion you do yourself.

`StripLoads::induced_angle` and `induced_velocity` are declared but never
filled. Per-panel induced velocity is in `SolutionState::induced_velocity`.

## Dimensional loads

`AerodynamicResults` is entirely dimensionless. Multiply by
`dynamic_pressure * reference_area`, and by `reference_chord` as well for
moments — including roll and yaw. For dimensional loads computed for you, with
right-hand moment signs, see [MeshTrimSolver](./mesh-trim-solver.md).

[`ForceCalculator` rustdoc](https://docs.rs/zephyr-vlm/latest/zephyr_vlm/aero/struct.ForceCalculator.html)
