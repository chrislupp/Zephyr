---
sidebar_position: 3
---

# Control Surfaces

A deflection in Zephyr rotates the **flow boundary condition**, not the mesh.
Panels affected by a command have their normals rotated about the hinge axis;
the lattice geometry and the wake stay where they were. That is the standard
quasi-steady VLM treatment, it is linear in nothing (the rotation is exact), and
it is what the control derivatives differentiate.

Deflections live in `FlightCondition::control_deflections`, an ordered list of
`(name, radians)` pairs:

```rust
use zephyr_vlm::solver::FlightCondition;

let condition = FlightCondition {
    alpha: 4.0_f64.to_radians(),
    control_deflections: vec![
        ("elevator".into(), -3.0_f64.to_radians()),
        ("aileron".into(), 5.0_f64.to_radians()),
    ],
    ..FlightCondition::default()
};
```

It is a list rather than a map on purpose: rotations compose, so the order
matters. Commands are applied in list order, and within one panel in definition
order. Keep that order fixed across an optimization or a sweep, or you are
differentiating a different function each iteration.

An omitted command is zero. An unknown name is an error, not a no-op.

## Section geometry: `ControlSurfaceDefinition`

This is AVL's `CONTROL`, attached to the sections that bound the deflecting
span.

```rust
use zephyr_vlm::core::Vec3;
use zephyr_vlm::geometry::{ControlSurfaceDefinition, Section};

let mut inboard = Section::new(Vec3::new(4.0, 1.0, 0.0), 0.6, 0.0);
let mut outboard = Section::new(Vec3::new(4.0, 3.0, 0.0), 0.5, 0.0);

let elevator = ControlSurfaceDefinition {
    name: "elevator".into(),
    hinge_x_fraction: 0.7,
    hinge_vector: Vec3::zeros(),
    sign: 1.0,
    gain: 1.0,
};
inboard.add_control_surface(elevator.clone());
outboard.add_control_surface(elevator);
```

| Field | Meaning |
| --- | --- |
| `name` | the command. Definitions sharing a name across surfaces are one command. |
| `hinge_x_fraction` | hinge position as a fraction of local chord. **Negative** means a leading-edge control: panels *ahead* of the hinge deflect. |
| `hinge_vector` | hinge axis. Zero falls back to the strip's spanwise unit vector, which is what you want for an ordinary trailing-edge flap. |
| `sign` | the deflection sign on the mirrored copy: `+1` symmetric (flap, elevator), `-1` antisymmetric (aileron, rudder on a duplicated fin). |
| `gain` | rotation per radian of command. |

Two things about this interface routinely surprise people:

**Both bounding sections need the definition.** A command is realized on the
interval between two sections only if *both* carry a definition with that name.
A definition on one section alone is silently inert for that interval. Hinge
fraction and gain are then interpolated across the interval.

**Partial panels get partial gain.** A panel straddling the hinge line receives
a gain scaled by the fraction of its chordwise extent that lies aft of the
hinge. This keeps the response continuous as the hinge moves, but it also means
the effective hinge is smeared over one panel. Align a chordwise panel edge with
the hinge when the derivative matters.

## Direct meshes: `MeshControlSurface`

For a direct mesh there is no chord fraction to speak of, because the mesh may
already be deformed. Hinges are therefore given in **material mesh
coordinates**.

```rust
use zephyr_vlm::geometry::MeshControlSurface;

// A 4 x 6 panel mesh: strips 1 through 5, hinge at chordwise node row 3.
let flap = MeshControlSurface::trailing_edge("flap", 1..6, vec![3.0; 7]);
// wing.controls.push(flap);
```

| Field | Meaning |
| --- | --- |
| `name` | the command, shared across surfaces as above. |
| `span` | active spanwise strips, start inclusive, end exclusive. |
| `hinge` | one coordinate per spanwise **node** station, so `ny + 1` entries, including inactive ones. |
| `leading_edge` | rotate the panels ahead of the hinge instead of behind it. |
| `gain` | rotation per radian, right-hand rule about the hinge running from station `j` to `j + 1`. |
| `duplicate_sign` | `+1` symmetric, `-1` antisymmetric on the mirrored copy. |

A hinge coordinate is an index into the chordwise node rows: `0` is the leading
edge, `nx` the trailing edge, and `2.5` is halfway between rows 2 and 3. It is
**material**, not geometric — it names a position in the grid, not a fraction of
the deformed chord. The consequence is the useful one: the hinge rides along
with bending and twist, and no re-derivation of a chord fraction is needed when
the structure deflects.

Coverage of a partial panel uses the mean overlap at the panel's two span
stations, the same smearing as the section interface.

Mesh derivatives hold the material hinge coordinates and the active span
interval fixed. Hinge placement and span extent are configuration, not design
variables; there are no derivatives with respect to them.

## Mixing and linking commands

- **One command, several surfaces.** Give matching `name`s. An elevator split
  across two tail panels is one command.
- **Independent left and right.** Use distinct names on separate meshes or
  surfaces. This is how you get an aileron pair you can command asymmetrically,
  rather than relying on `duplicate_sign`.
- **Several commands on one panel.** Add several definitions. A panel can be
  both flap and aileron; the rotations compose.
- **Fixed linear mixers.** For trim, a `trim::ControlMixer` reduces redundant
  physical surfaces to independent commands through
  `delta = offsets + matrix * commands`. See [Trim](./trim.md).

Repeated entries of the same name in one `FlightCondition` list are legal for
backward compatibility. A named derivative then varies every entry of that name
together, and the derivative of an omitted command appends it at zero.

## Deflected or deformed — pick one

This matters when an external code owns the mesh. Hand Zephyr the
structurally deformed mesh **before** control actuation and put the deflection
in the flight condition. Zephyr rotates the normals.

If your structural model has already physically rotated the control-surface
nodes, do **not** also command the deflection here — you would apply it twice.
Supply the deflected mesh and omit the command.

## Derivatives

Control derivatives are analytical, from the same adjoint machinery as
everything else. Nothing is finite-differenced.

| Call | Returns |
| --- | --- |
| `GradientCalculator::compute_control_gradients` | coefficient and far-field induced-drag derivatives per radian, for every named command, honoring `OutputSelection` |
| `StabilityDerivativeCalculator::compute` | per-command `ControlDerivatives` alongside the angle and rate derivatives |
| `StripLoadAdjoint::control_cotangents` | derivatives of an arbitrary functional of the distributed strip forces and moments |
| `analysis::Input::Control(name)` | one requested `(output, command)` pair, in a selective request |

Mesh-node and section-geometry gradients also work at nonzero deflection, and
the mesh gradients include the hinge axis moving with the mesh.

Hinge moments are not computed: `ControlDerivatives::hinge_moment` is present in
the struct and is reported as zero.

The native example is `cargo run --example control_surfaces`; see
[Control Surfaces (example)](../examples/control-surfaces.md).
