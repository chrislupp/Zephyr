---
sidebar_position: 3
---

# Surface

`geometry::Surface` is one lifting component: a list of spanwise `Section`s plus
the discretization, transforms, duplication, flags, and drag polar that apply to
all of them.

A surface needs at least two sections and a nonzero projected span, or the
mesher produces nothing for it.

## Discretization

| Accessor | Meaning |
| --- | --- |
| `chordwise_panel_count`, `spanwise_panel_count` | panel counts |
| `chordwise_spacing`, `spanwise_spacing` | `SpacingType`; `Cosine` by default in both directions |
| `set_numeric_chordwise_spacing`, `set_numeric_spanwise_spacing` | AVL's signed convention: zero uniform, positive cosine, negative trailing-edge sine |
| `spanwise_strip_count()` | the strip count the mesher will actually build |

The spanwise total is divided among section intervals in proportion to projected
arc length, truncated, with a floor of one panel per interval, so the realized
total can be below the requested one. `spanwise_strip_count()` is computed from
the mesher's own layout and is the number to trust.

## Placement

| Accessor | Meaning |
| --- | --- |
| `translation` | per-axis offset, AVL's `TRANSLATE` |
| `scale` | per-axis scale factors, AVL's `SCALE` |
| `angle_offset_degrees` | incidence added to every section, AVL's `ANGLE`, in **degrees** |
| `y_duplicate`, `set_y_duplicate`, `clear_y_duplicate` | mirror plane, AVL's `YDUPLICATE` |
| `component_index` | component grouping for core regularization |

Duplication generates a real mirrored surface with its own panels and unknowns,
placed immediately after the original in the lattice. It is the duplication the
adjoints and the symmetry accelerations require.

Component grouping matters for junctions: the regularization core is smaller
between panels of the same component than across components, which is what keeps
a deliberately touching wing-and-flap pair from going singular.

## Flags

| Accessor | Status |
| --- | --- |
| `affected_by_alpha_beta` | **active**; false gives this surface a zero right-hand side (AVL `NOALBE`) |
| `has_wake` | recorded per panel, not acted on by the solver (AVL `NOWAKE`) |
| `contributes_to_loads` | recorded, not applied by the force integration (AVL `NOLOAD`) |

## Drag polar

`cl_cd_polar` and `set_cl_cd_polar` hold AVL's `CDCL`: six values, three
$(C_l, C_d)$ pairs, as a surface-wide fallback for sections that carry no polar
of their own. `has_viscous_drag()` reports whether any is set.

The polar is a Lagrange quadratic through the three points and is deliberately
unclamped, so it extrapolates parabolically outside the fitted range.

A mesh carrying any section polar is **rejected** by the geometry and strip-load
adjoints; flight-condition and control derivatives carry the viscous chain
normally.

## Geometric reports

`planform_area()`, `span()`, and `mean_aerodynamic_chord()` are computed from
the sections. They are reports, not reference quantities — set the references on
the `Aircraft` yourself.

[`Surface` rustdoc](https://docs.rs/zephyr-vlm/latest/zephyr_vlm/geometry/struct.Surface.html)
