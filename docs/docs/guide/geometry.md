---
sidebar_position: 2
---

# Geometry

There are two ways to describe geometry, and they are not two spellings of the
same thing:

- **Sections.** An `Aircraft` of `Surface`s, each a list of spanwise `Section`s
  with chords, twists, and airfoils. Zephyr meshes it. This is the AVL-style
  description, and it is what the AVL reader produces.
- **Direct meshes.** A `MeshInput` of `SurfaceMesh` node grids you supply
  yourself. Zephyr panels exactly what you hand it. This is the interface for
  aeroelasticity and for optimizers that own their own parameterization.

Both feed the same solver. Pick one per model: `VortexLatticeSolver::set_aircraft`
and `set_mesh` each replace whatever the other established.

## Surfaces and sections

A `Surface` needs at least two sections and a nonzero projected span. Each
`Section` carries a leading-edge point, a chord, a twist in **degrees**, and
optionally an airfoil, a lift-slope factor, a viscous polar, spanwise overrides,
and control-surface and design-variable metadata.

```rust
use zephyr_vlm::core::{SpacingType, Vec3};
use zephyr_vlm::geometry::{Airfoil, Section, Surface};

let mut wing = Surface::new("Wing");
wing.set_chordwise_panel_count(8);
wing.set_spanwise_panel_count(24);
wing.set_chordwise_spacing(SpacingType::Cosine);
wing.set_spanwise_spacing(SpacingType::Uniform);
wing.set_y_duplicate(0.0);

let mut root = Section::new(Vec3::new(0.0, 0.0, 0.0), 2.0, 2.0);
root.set_airfoil(Airfoil::naca4(2412, 61)?);
wing.add_section(root);

let mut tip = Section::new(Vec3::new(1.8, 7.0, 0.3), 0.8, -2.0);
tip.set_airfoil(Airfoil::naca4(2410, 61)?);
wing.add_section(tip);
```

Sections are added in increasing span order. The leading-edge points define
sweep and dihedral. Between two sections, the leading-edge position and the
chord are interpolated linearly in projected arc length. Twist follows AVL's
chord-weighted convention rather than a linear blend of angles: the strip twist
is the arctangent of the interpolated $c\sin\theta$ over the interpolated
$c\cos\theta$. Camber slope, viscous polar coefficients, and the lift-slope
factor are blended across the interval as well.

Panels are flat in the vertical sense: a strip's bound vortices and control
points are placed along the chord line at the AVL chordwise stations, and camber
enters through the panel normal, not by displacing the lattice onto the camber
surface.

## Panel spacing

`SpacingType` chooses the distribution of panel edges:

| Variant | Distribution |
| --- | --- |
| `Uniform` | equal intervals |
| `Cosine` | clustered at both ends |
| `Sine` | clustered at the trailing end |
| `Blend` | the average of `Uniform` and `Cosine` |
| `None` | no override — inherit from the surface |

`Surface::default()` is `Cosine` in both directions. A `Section` may override
the spanwise distribution for the interval that begins at it; `None` there falls
back to the surface's own setting.

Cosine spanwise spacing is usually right for an isolated wing, where the
circulation gradient is steepest at the tip. Uniform spacing is usually right
when a strip boundary has to line up with something — a control-surface edge, a
structural rib, or another surface's panel edges.

AVL's signed numeric spacing parameters are accepted directly through
`set_numeric_spanwise_spacing` and `set_numeric_chordwise_spacing`: zero means
uniform, positive means cosine, negative means trailing-edge sine. The AVL
reader uses these, so a parsed file keeps its spacing.

### How the spanwise count is distributed

Set `Section::set_spanwise_panel_count` and that interval gets exactly that many
panels. Leave it at zero and the surface's total is divided among the intervals
in proportion to their projected arc length, truncated toward zero, with a
floor of one panel per interval. A surface with no count at all gets ten.

Truncation means the realized total can be lower than the number you asked for.
`Surface::spanwise_strip_count` reports what the mesher will actually build, and
it is computed from the same layout the mesher uses — so use it rather than
recomputing a share yourself.

## Duplication

`set_y_duplicate(plane)` generates a mirrored copy of the surface about that `y`
plane. The copy is real: it has its own panels, its own unknowns, and its own
entry in `VortexLattice::surface_info`, immediately after the original.
`SurfaceMeshInfo::original_surface_index` points back at the surface it was
mirrored from.

This is the duplication to use. The alternative — image planes through
`VortexLatticeSolver::set_symmetry` — adds reflected influence without adding
unknowns, which is cheaper but is not supported by the geometry or strip-load
adjoints. If you intend to differentiate the geometry, duplicate explicitly.

A duplicated half-model is also the precondition for the three symmetry
accelerations: mirror reuse, symmetric factorization, and half-wing mode. See
[Performance](./performance.md).

## Airfoils

`Airfoil` holds a camber line and a thickness distribution. Only the camber
matters to the solve: the panel normal is rotated by the local camber slope.
Thickness is kept for visualization and for reporting.

| Constructor | Source |
| --- | --- |
| `Airfoil::flat_plate()` | zero camber, zero thickness |
| `Airfoil::naca4(2412, 61)` | four-digit NACA from a numeric code |
| `Airfoil::naca4_str("2412", 61)` | the same from text |
| `Airfoil::from_coordinate_file(path)` | Selig or Lednicer coordinate file |
| `Airfoil::from_coordinates(&x, &y, is_camber_line)` | explicit camber line, or closed coordinates to split |

A section with no airfoil behaves as a flat plate. `Section::airfoil_mut`
creates one on demand.

`Section::set_cl_alpha_factor` is AVL's `CLAF`: it scales the section lift slope
by moving the chordwise control-point position, so a value above one raises
$C_{L\alpha}$ for that section. It is applied in the chordwise spacing
construction, so it is exact for every spacing law.

## Surface flags

Three per-surface flags correspond to AVL's `NOALBE`, `NOWAKE`, and `NOLOAD`.
Only the first currently changes a result:

| Setter | Status |
| --- | --- |
| `set_affected_by_alpha_beta(false)` | **Active.** The freestream does not drive this surface's boundary condition: its right-hand side is zero. |
| `set_has_wake(false)` | Recorded on each generated panel as `enforce_flow_tangency = false`. The solver does not currently act on it. |
| `set_contributes_to_loads(false)` | Recorded on the surface. The force integration does not currently exclude it. |

A `NOALBE` surface still induces velocity on everything else, still carries
loads, and can still be deflected by a control; it simply does not respond to
the freestream itself. The usual use is a fuselage or nacelle carried as a
lifting shell.

Because the other two flags are inert, a surface you want excluded from the
reported loads has to be excluded by you — the per-surface records in
`AerodynamicResults::surface_loads` make that a subtraction rather than a
re-run. The AVL reader sets all three faithfully, so a parsed file keeps its
intent even where the solver does not yet honor it.

`set_component_index` groups surfaces that are physically one body. The core
regularization radius is smaller between panels of the same component than
across components, which is what keeps a wing-and-flap pair or an
intentionally-touching junction from producing a singular influence.

`set_translation` and `set_scale` apply a per-axis offset and scale to the whole
surface after its sections are read, matching AVL's `TRANSLATE` and `SCALE`.
`set_angle_offset_degrees` adds an incidence offset to every section, matching
`ANGLE`.

## Viscous drag

Two independent mechanisms, both additive onto the inviscid result:

- `Aircraft::set_profile_drag_coefficient` is a single constant added to the
  aircraft's drag coefficient. It is AVL's header `CDp`.
- `Surface::set_cl_cd_polar` and `Section::set_viscous_polar` are AVL's `CDCL`:
  three $(C_l, C_d)$ pairs, packed as six values. Each strip evaluates its own
  $C_d$ from its own local $C_l$ and contributes it to the drag.

The polar is a Lagrange quadratic through the three points, deliberately
unclamped, so it extrapolates parabolically outside the fitted range rather than
saturating. That is a choice you should know about before relying on it far from
where you fitted it. Strip viscous drag is applied along the freestream
direction, so it rotates with $\alpha$, and it is reported separately as
`AerodynamicResults::viscous_drag` as well as being folded into `drag`.

## Bodies

`Body` and `BodySection` describe elliptic-section bodies, and the AVL reader
populates them from `BODY` and `BSEC` blocks so that a parsed file round-trips
faithfully. **They do not contribute to the solution.** There is no source-line
or doublet-line body model: bodies are stored geometry, available to your own
code through `Aircraft::bodies`, and nothing more. Model a fuselage you need
aerodynamically as a lifting surface with `NOALBE`, or accept its absence.

## Direct meshes

A `SurfaceMesh` is an `(nx+1) × (ny+1)` grid of nodes, row-major, indexed
`i * (ny + 1) + j`, with `i` chordwise from leading to trailing edge and `j`
spanwise.

```rust
use zephyr_vlm::core::Vec3;
use zephyr_vlm::geometry::{MeshInput, SurfaceMesh};

let mut wing = SurfaceMesh::new("Wing", 4, 16);
wing.y_duplicated = true;
for i in 0..=wing.nx {
    for j in 0..=wing.ny {
        let x = i as f64 / wing.nx as f64;
        let y = 5.0 * j as f64 / wing.ny as f64;
        wing.set_node(i, j, Vec3::new(x, y, 0.0))?;
    }
}

let mut input = MeshInput::default();
input.surfaces.push(wing);
input.reference.area = 10.0;
input.reference.chord = 1.0;
input.reference.span = 10.0;
input.reference.moment_reference_point = Vec3::new(0.25, 0.0, 0.0);
input.validate()?;
```

Panels are built from the nodes directly: the bound vortex sits a quarter of the
way from node row `i` to row `i + 1`, the control point three quarters of the
way, and the normal comes from the cell's diagonals. So camber and twist are
whatever the node positions say they are — there is no airfoil to blend and no
twist angle to apply.

That is the point of the interface. A structural solver hands over deformed node
coordinates and the aerodynamics follows them exactly. It also means a direct
mesh and a section model of nominally the same wing will **not** produce
identical lattices: the section mesher places nodes on a camber line at AVL's
chordwise stations, while a direct mesh places them wherever you put them.

`SurfaceMesh` fields mirror the surface flags: `component_index`,
`affected_by_freestream`, `has_wake`, `y_duplicated`, `y_duplicate_plane`, and a
`controls` list of [mesh control surfaces](./control-surfaces.md).
`validate()` checks the dimensions and that every coordinate is finite;
`MeshInput::validate` checks every surface and the reference quantities.

Reference quantities travel in `MeshInput::reference`, not on the surfaces.
`ReferenceQuantities::from(&aircraft)` extracts them from a section model, which
is how you move a model from one interface to the other without changing its
normalization.
