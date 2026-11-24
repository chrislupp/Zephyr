---
sidebar_position: 4
---

# VTK Visualization

Four VTK XML products written from one solution, with no system VTK library
involved.

```bash
cargo run --example vtk_visualization -- vtk_output
```

Source: `examples/vtk_visualization.rs`. The argument is the output directory
and defaults to `vtk_output`.

## What it writes

| File | Format | Contents |
| --- | --- | --- |
| `surface_mesh.vtu` | unstructured grid | the panel mesh with per-cell fields |
| `spanwise_loads.vtp` | polydata | the spanwise load distribution as lines |
| `wake.vtp` | polydata | trailing wake lines |
| `velocity_field.vtp` | polydata | induced velocity at the panel points |

`VtkExporter::export_all` writes all four into a directory it creates. The
individual `export_surface_mesh`, `export_spanwise_loads`, `export_wake`, and
`export_velocity_field` functions take a file path each.

Every export is written to a temporary sibling file and renamed only after
serialization succeeds, so a failed export never leaves a half-written file for
ParaView to choke on.

## What it does

```rust
use zephyr_vlm::aero::ForceCalculator;
use zephyr_vlm::core::{DEGREES_TO_RADIANS, Vec3};
use zephyr_vlm::geometry::{Aircraft, Section, Surface};
use zephyr_vlm::solver::{FlightCondition, VortexLatticeSolver};
use zephyr_vlm::visualization::{VtkExportOptions, VtkExporter, is_available};

fn main() -> zephyr_vlm::Result<()> {
    if !is_available() {
        println!("VTK support is disabled in this build");
        return Ok(());
    }

    let mut wing = Surface::new("Tapered wing");
    wing.set_chordwise_panel_count(6);
    wing.set_spanwise_panel_count(16);
    wing.set_y_duplicate(0.0);
    wing.add_section(Section::new(Vec3::new(0.0, 0.0, 0.0), 1.0, 2.0));
    wing.add_section(Section::new(Vec3::new(0.25, 5.0, 0.0), 0.5, -1.0));

    let mut aircraft = Aircraft::new("VTK example");
    aircraft.set_reference_area(7.5);
    aircraft.set_reference_chord(0.75);
    aircraft.set_reference_span(10.0);
    aircraft.add_surface(wing);

    let condition = FlightCondition {
        alpha: 5.0 * DEGREES_TO_RADIANS,
        ..FlightCondition::default()
    };
    let mut solver = VortexLatticeSolver::new();
    solver.set_aircraft(&aircraft)?;
    let solution = solver.solve(&condition)?;
    let results = ForceCalculator::new().compute_for_aircraft(
        &aircraft, solver.mesh(), &solution, &condition,
    );

    VtkExporter::export_all(
        "vtk_output",
        solver.mesh(),
        &solution,
        &results,
        &condition,
        aircraft.reference_chord(),
        VtkExportOptions::default(),
    )?;
    Ok(())
}
```

The wing is a half model with `set_y_duplicate(0.0)`, root at `y = 0`, tapered
and twisted from $+2^\circ$ at the root to $-1^\circ$ at the tip. The export
includes the mirrored half, because the lattice does.

## Feature gate

VTK export is behind the default `vtk` feature. In a build without it,
`is_available()` returns `false` and every export function returns
`ZephyrError::FeatureDisabled` — which is why the example checks first instead of
failing. The underlying crate is the pure-Rust `vtkio`; there is no system VTK
dependency either way.

## Options

```rust
let options = VtkExportOptions {
    export_normals: true,
    export_indices: true,
    wake_length_chords: 20.0,
    binary_format: false,
    ..VtkExportOptions::default()
};
```

| Field | Default | Effect |
| --- | --- | --- |
| `export_cp` | `true` | per-panel pressure difference |
| `export_circulation` | `true` | per-panel circulation |
| `export_normals` | `false` | panel normal vectors |
| `export_indices` | `false` | the `Surface` and `Strip` index arrays — useful for debugging ordering |
| `export_wake` | `true` | include the wake product |
| `wake_length_chords` | `10.0` | wake length in reference chords |
| `export_induced_velocity` | `true` | include the velocity-field product |
| `export_spanwise_distributions` | `true` | include the spanwise-loads product |
| `binary_format` | `true` | base64 binary arrays rather than ASCII |

Binary is the default and is what you want for anything large. ASCII is
readable, diffable, and worth switching to when you are checking a field by eye.
Both parse back with `vtkio`.

Two arrays in the spanwise-loads product are currently always zero —
`InducedAlpha` and its `InducedVelocity`, because `StripLoads::induced_angle`
and `StripLoads::induced_velocity` are never filled. Per-panel induced velocity
*is* real, and reaches the velocity-field product's `InducedVelocity` from
`SolutionState::induced_velocity`.

Array names, for colouring in a viewer: the surface mesh carries `Cp`,
`Circulation`, `Normal`, `Surface`, and `Strip` as cell data; the spanwise loads
carry `CL`, `CD`, `CY`, `SpanLoading`, `CM`, and `Surface` as point data; the
velocity field carries `InducedVelocity`, `TotalVelocity`, and
`VelocityMagnitude`.

## Viewing

The files open in ParaView or VisIt. The surface mesh is the one to colour by
`Cp` or `Circulation`; the spanwise loads are lines you plot along the span; the
wake shows the trailing geometry the Trefftz integration sees.
