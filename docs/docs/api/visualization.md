---
sidebar_position: 16
---

# VtkExporter

`visualization::VtkExporter` writes VTK XML products through the pure-Rust
`vtkio` crate. There is no system VTK dependency.

Behind the default `vtk` feature. In a build without it, `is_available()`
returns `false` and every export returns
`ZephyrError::FeatureDisabled("vtk")` — so check first rather than failing.

The worked example is [VTK Visualization](../examples/vtk-visualization.md).

```rust
use zephyr_vlm::visualization::{VtkExportOptions, VtkExporter, is_available};

if is_available() {
    VtkExporter::export_all(
        "vtk_output",
        solver.mesh(),
        solution,
        results,
        condition,
        1.0, // reference chord, for the wake length
        VtkExportOptions::default(),
    )?;
}
```

## Products

| Method | File | Contents |
| --- | --- | --- |
| `export_surface_mesh` | `.vtu` | the panel mesh, cell data |
| `export_spanwise_loads` | `.vtp` | the spanwise load distribution as lines, point data |
| `export_wake` | `.vtp` | trailing wake lines |
| `export_velocity_field` | `.vtp` | induced and total velocity at the panel points |
| `export_all` | a directory | all four, as `surface_mesh.vtu`, `spanwise_loads.vtp`, `wake.vtp`, `velocity_field.vtp` |

`export_all` creates the directory. The individual functions each take a file
path.

Every export writes to a temporary sibling and renames only after serialization
succeeds, so a failed export never leaves a half-written file behind.

An invalid solution or an empty lattice is `ZephyrError::InvalidSolution`.

## Array names

| Product | Arrays |
| --- | --- |
| Surface mesh | `Cp`, `Circulation`, `Normal`, `Surface`, `Strip` (cell data) |
| Spanwise loads | `CL`, `CD`, `CY`, `SpanLoading`, `CM`, `InducedAlpha`, `InducedVelocity`, `Surface` (point data) |
| Velocity field | `InducedVelocity`, `TotalVelocity`, `VelocityMagnitude`, `Surface` |

The spanwise product's `InducedAlpha` and `InducedVelocity` are currently
**always zero**, because `StripLoads::induced_angle` and
`StripLoads::induced_velocity` are never filled. The velocity field's
`InducedVelocity` is real, from `SolutionState::induced_velocity`.

## `VtkExportOptions`

| Field | Default | Effect |
| --- | --- | --- |
| `export_cp` | `true` | per-panel pressure difference |
| `export_circulation` | `true` | per-panel circulation |
| `export_normals` | `false` | panel normal vectors |
| `export_indices` | `false` | the `Surface` and `Strip` index arrays |
| `export_wake` | `true` | include the wake product |
| `wake_length_chords` | `10.0` | wake length in reference chords |
| `export_induced_velocity` | `true` | include the velocity-field product |
| `export_spanwise_distributions` | `true` | include the spanwise-loads product |
| `binary_format` | `true` | base64 binary arrays rather than ASCII |

Binary is the default. ASCII is readable and diffable, and worth switching to
when checking a field by eye; both parse back with `vtkio`.

The index arrays are the quickest way to confirm the panel and strip ordering
described in [Conventions](../getting-started/conventions.md).

[`VtkExporter` rustdoc](https://docs.rs/zephyr-vlm/latest/zephyr_vlm/visualization/struct.VtkExporter.html)
