---
sidebar_position: 4
---

# Section

`geometry::Section` is one spanwise station: a leading-edge point, a chord, a
twist, and optionally an airfoil, a lift-slope factor, a viscous polar,
spanwise overrides, and control and design metadata.

```rust
use zephyr_vlm::core::Vec3;
use zephyr_vlm::geometry::{Airfoil, Section};

let mut root = Section::new(Vec3::new(0.0, 0.0, 0.0), 2.0, 2.0);
root.set_airfoil(Airfoil::naca4(2412, 61)?);
```

`Section::new(leading_edge, chord, twist_degrees)` — **twist is in degrees**.
`twist_degrees()` and `twist_radians()` both read it; only the degree setter
exists, and section twist derivatives are reported per degree.

`Section::with_airfoil` constructs one with an airfoil in a single call.

## Interpolation between stations

Between two sections, the leading-edge position and the chord are interpolated
linearly in projected arc length. Twist follows AVL's chord-weighted convention
rather than a linear blend of angles: the strip twist is the arctangent of the
interpolated $c\sin\theta$ over the interpolated $c\cos\theta$. Camber slope,
viscous polar coefficients, and the lift-slope factor are blended across the
interval too.

Sections are added in increasing span order.

## Airfoil

| Accessor | Meaning |
| --- | --- |
| `airfoil()` | the airfoil, or a shared flat plate when none was set |
| `airfoil_mut()` | the airfoil, creating a flat plate if needed |
| `set_airfoil`, `has_airfoil` | replace, test |

Only the **camber** reaches the solve, as a rotation of the panel normal.
Thickness is carried for visualization and reporting.

`cl_alpha_factor` is AVL's `CLAF`: it scales the section lift slope by moving
the chordwise control-point position, applied inside the chordwise spacing
construction so it is exact for every spacing law.

## Overrides

| Accessor | Meaning |
| --- | --- |
| `spanwise_panel_count` | panels in the interval beginning at this section; zero means inherit |
| `spanwise_spacing` | distribution for that interval; `None` means inherit from the surface |
| `set_numeric_spanwise_spacing` | AVL's signed convention |
| `viscous_polar` | this section's `CDCL`, six values |

## Metadata

`add_control_surface` and `control_surfaces()` hold
`ControlSurfaceDefinition`s. A command is realized on an interval only when
**both** bounding sections carry a definition with that name — a definition on
one section alone is silently inert. Fields:

| Field | Meaning |
| --- | --- |
| `name` | the command; matching names across surfaces are one command |
| `hinge_x_fraction` | hinge as a chord fraction; **negative** means a leading-edge control |
| `hinge_vector` | hinge axis; zero falls back to the strip's spanwise unit vector |
| `sign` | deflection sign on the mirrored copy: `+1` symmetric, `-1` antisymmetric |
| `gain` | rotation per radian of command |

A panel straddling the hinge gets a gain scaled by its chordwise overlap, so the
effective hinge is smeared over one panel. Align a panel edge with the hinge
when the derivative matters.

`add_design_variable` and `design_variables()` hold `DesignVariableDefinition`
— a name and a station weight, parsed from AVL's `DESIGN`. This is metadata the
library records for you; the gradient interface works from
[`GeometryVariable`](./gradient-calculator.md) rather than from these names.

[`Section` rustdoc](https://docs.rs/zephyr-vlm/latest/zephyr_vlm/geometry/struct.Section.html)
