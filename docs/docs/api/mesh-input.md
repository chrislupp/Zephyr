---
sidebar_position: 5
---

# MeshInput and SurfaceMesh

The direct-mesh interface. `geometry::MeshInput` is a list of `SurfaceMesh` node
grids plus the reference quantities; the solver panels exactly what it is given.

This is the interface for aeroelasticity and for optimizers that own their own
parameterization. See [Coupled Interfaces](../guide/coupled-interfaces.md).

## `SurfaceMesh`

An `(nx+1) × (ny+1)` grid of nodes, **row major**, indexed `i * (ny + 1) + j`,
with `i` chordwise from leading to trailing edge and `j` spanwise.

```rust
use zephyr_vlm::core::Vec3;
use zephyr_vlm::geometry::SurfaceMesh;

let mut wing = SurfaceMesh::new("Wing", 4, 16);
wing.y_duplicated = true;
for i in 0..=wing.nx {
    for j in 0..=wing.ny {
        wing.set_node(i, j, Vec3::new(i as f64 / 4.0, 5.0 * j as f64 / 16.0, 0.0))?;
    }
}
wing.validate()?;
```

| Field | Meaning |
| --- | --- |
| `name` | surface name, carried into the lattice and the load records |
| `nx`, `ny` | chordwise and spanwise **panel** counts |
| `nodes` | the row-major node list |
| `controls` | `MeshControlSurface` definitions, applied before actuation |
| `component_index` | component grouping for core regularization |
| `affected_by_freestream` | false gives a zero right-hand side |
| `has_wake` | recorded; the solver does not currently act on it |
| `y_duplicated`, `y_duplicate_plane` | generate a mirrored copy |

`node(i, j)` and `set_node(i, j, point)` are the checked accessors;
`node_count()` and `panel_count()` are const. `validate()` checks the dimensions
and that every coordinate is finite.

## How a mesh becomes panels

Panels come straight from the nodes. The bound vortex sits a quarter of the way
from node row `i` to row `i + 1`, the control point three quarters of the way,
and the panel normal comes from the cell's diagonals.

So **camber and twist are the node positions**. There is no airfoil to blend and
no twist angle to apply. That is what makes a deformed mesh from a structural
solver mean what the structural solver intended — and it is also why a direct
mesh and a section model of nominally the same wing do not produce identical
lattices.

## `ReferenceQuantities`

```rust
use zephyr_vlm::core::Vec3;
use zephyr_vlm::geometry::{MeshInput, ReferenceQuantities};

let mut input = MeshInput::default();
input.reference = ReferenceQuantities {
    area: 10.0,
    chord: 1.0,
    span: 10.0,
    moment_reference_point: Vec3::new(0.25, 0.0, 0.0),
};
input.validate()?;
```

References live on the `MeshInput`, not on the surfaces.
`ReferenceQuantities::from(&aircraft)` extracts them from a section model, which
is how a model moves between the two interfaces without changing its
normalization.

`MeshInput::validate()` validates every surface and the reference values.
`set_mesh` validates for you and reports `ZephyrError::InvalidMesh`.

## `MeshControlSurface`

Hinges in **material mesh coordinates**: `0` is the leading edge, `nx` the
trailing edge, and `2.5` is halfway between node rows 2 and 3. The coordinate
names a position in the grid, not a fraction of the deformed chord, so the hinge
rides along with bending and twist.

```rust
use zephyr_vlm::geometry::MeshControlSurface;

// 4 x 6 panel mesh: strips 1 through 5, hinge at chordwise node row 3.
let flap = MeshControlSurface::trailing_edge("flap", 1..6, vec![3.0; 7]);
```

| Field | Meaning |
| --- | --- |
| `name` | the command; matching names across surfaces are one command |
| `span` | active strips, start inclusive, end exclusive |
| `hinge` | one coordinate per spanwise **node** station: `ny + 1` entries, inactive ones included |
| `leading_edge` | rotate panels ahead of the hinge instead of behind it |
| `gain` | rotation per radian, right-hand rule about the hinge from station `j` to `j + 1` |
| `duplicate_sign` | `+1` symmetric, `-1` antisymmetric on the mirrored copy |

Mesh derivatives hold the material hinge coordinates and the active span
interval fixed; there are no derivatives with respect to them. Partial coverage
uses the mean overlap at the panel's two span stations.

See [Control Surfaces](../guide/control-surfaces.md) for the deflected-or-deformed
rule, which is the one that produces wrong answers quietly if you get it
backwards.

[`MeshInput` rustdoc](https://docs.rs/zephyr-vlm/latest/zephyr_vlm/geometry/struct.MeshInput.html)
