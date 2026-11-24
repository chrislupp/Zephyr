---
sidebar_position: 2
---

# Aircraft

`geometry::Aircraft` is the root of a section model: surfaces, bodies, and the
reference quantities that normalize every reported coefficient.

```rust
use zephyr_vlm::core::Vec3;
use zephyr_vlm::geometry::{Aircraft, Surface};

let wing = Surface::new("Wing"); // sections added elsewhere

let mut aircraft = Aircraft::new("Example");
aircraft.set_reference_area(10.0);
aircraft.set_reference_chord(1.0);
aircraft.set_reference_span(10.0);
aircraft.set_moment_reference_point(Vec3::new(0.25, 0.0, 0.0));
aircraft.add_surface(wing);
```

## Reference quantities

| Accessor | Used for |
| --- | --- |
| `reference_area` | every force and moment coefficient |
| `reference_chord` | **every** moment coefficient, roll and yaw included |
| `reference_span` | aspect ratio and span efficiency only |
| `moment_reference_point` | the point all moments are taken about |

Nothing is inferred from the planform. A wrong reference area silently rescales
every coefficient and is the most common first mistake;
`total_planform_area()` and `total_span()` are available to check against.

## Flow and symmetry

| Accessor | Meaning |
| --- | --- |
| `mach_number` | adopted as the solver's current Mach by `set_aircraft` |
| `profile_drag_coefficient` | AVL's header `CDp`, one constant added to drag |
| `y_symmetry`, `z_symmetry`, `z_symmetry_plane` | image-plane symmetry, as AVL's `iYsym`/`iZsym`/`Zsym` |

Image-plane symmetry is **not supported** by the geometry or strip-load
adjoints. Use `Surface::set_y_duplicate` for anything you intend to
differentiate — see [Geometry](../guide/geometry.md).

## Collections

| Accessor | Returns |
| --- | --- |
| `surfaces()`, `surfaces_mut()` | the lifting surfaces |
| `surface_named(name)` | a surface by name, checked |
| `bodies()` | body geometry, which **does not affect the solution** |
| `control_surface_names()` | ordered unique control names |
| `design_variable_names()` | ordered unique design-variable names |

## Aggregates

`total_vortex_count()`, `total_strip_count()`, `total_planform_area()`, and
`total_span()` report what the mesher will actually build, computed from the
same spanwise layout the mesher uses rather than from an independent count.

Mutating a surface through `surfaces_mut()` after a solver has been given this
aircraft does not invalidate that solver's caches. Call `set_aircraft` again, or
use `VortexLatticeSolver::mutate_mesh` for an in-place lattice change.

[`Aircraft` rustdoc](https://docs.rs/zephyr-vlm/latest/zephyr_vlm/geometry/struct.Aircraft.html)
