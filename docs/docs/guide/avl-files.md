---
sidebar_position: 9
---

# AVL Input Files

Zephyr reads AVL geometry (`.avl`) and mass (`.mass`) files, so an existing
model can be analyzed without being transcribed. The readers produce ordinary
`Aircraft` and `MassProperties` values — there is no separate "AVL mode", and
nothing downstream knows where the geometry came from.

```rust
use zephyr_vlm::io::{AvlReader, MassReader};

let aircraft = AvlReader::new().read_path("aircraft.avl")?;
let mass = MassReader::read_path("aircraft.mass")?.total_properties();
```

| Entry point | Reads |
| --- | --- |
| `AvlReader::read_path` | a file, resolving relative `AFILE` paths against its directory |
| `AvlReader::read_str` | text, with no base directory for airfoil files |
| `AvlReader::read_with_base` | text, with an explicit base directory |
| `MassReader::read_path` | a mass file |
| `MassReader::read_str` | mass-file text |

`AvlReader::set_verbose(true)` prints the configuration title and each surface
name to standard error as it reads.

## Geometry file

Blank lines and lines beginning with `#` or `!` are dropped before parsing.
Keywords are matched on their first four characters, case-insensitively, as AVL
does, so `SURFACE`, `SURF`, and `Surface` are the same token.

The header is positional, in this order:

1. configuration title
2. Mach number
3. three symmetry values: `iYsym`, `iZsym`, `Zsym`
4. three reference values: `Sref`, `Cref`, `Bref`
5. three moment reference coordinates: `Xref`, `Yref`, `Zref`
6. optionally `CDp`, a constant profile-drag coefficient — recognized only if
   the next line is not a `SURFACE` or `BODY` keyword

Then any number of `SURFACE` and `BODY` blocks.

### Surface keywords

A `SURFACE` block opens with its name on the following line and then a paneling
line: `Nchord Cspace [Nspan Sspace]`.

| Keyword | Effect |
| --- | --- |
| `SCALE` | per-axis scale factors |
| `TRANSLATE` | per-axis offset |
| `ANGLE` | incidence offset added to every section, in degrees |
| `YDUPLICATE` | mirror the surface about the given `y` plane |
| `COMPONENT` | component grouping for core regularization. AVL's `INDEX` synonym is **not** recognized; use `COMPONENT`. |
| `NOWAKE`, `NOLOAD`, `NOALBE` | the three surface flags |
| `CDCL` | six values: three $(C_l, C_d)$ polar points for the whole surface |
| `SECTION` | begin a spanwise station |

### Section keywords

A `SECTION` line carries `Xle Yle Zle chord incidence`, and optionally
`Nspan Sspace` as the sixth and seventh values to override the spanwise
discretization of the interval that begins at this section. Incidence is in
**degrees**.

| Keyword | Effect |
| --- | --- |
| `NACA` | four-digit code; the camber line is generated with 50 points |
| `AFILE` / `AIRFOIL` | coordinate file; the name runs to the end of the line, with any trailing comment stripped |
| `CLAF` | section lift-slope factor |
| `CDCL` | six values: this section's viscous polar |
| `CONTROL` | `name gain hinge_x hx hy hz [sign]`; `sign` defaults to 1 |
| `DESIGN` | `name weight`, recorded as design-variable metadata |

A control definition is realized only on intervals where **both** bounding
sections carry it — see [Control Surfaces](./control-surfaces.md).

### Spacing numbers

AVL's signed spacing parameters map as follows:

| Value | Distribution |
| --- | --- |
| magnitude below 0.01 | uniform |
| positive | cosine, clustered at both ends |
| negative | sine, clustered at the trailing end |

### Body blocks

`BODY` blocks are parsed — name, `SCALE`, `TRANSLATE`, `YDUPLICATE`, and `BSEC`
sections — and stored on the `Aircraft`. **They do not contribute to the
solution.** There is no body source or doublet model; see
[Geometry](./geometry.md).

### Deviations worth knowing

- A failed `AFILE` load is an error, not a silent fall back to a flat plate.
  Degrading camber while reporting success would change the flow the solve runs
  on.
- `NOWAKE` and `NOLOAD` are recorded but do not currently change a result.
- Mach number from the header becomes the solver's current Mach when the
  aircraft is installed with `set_aircraft`.
- Rolling and yawing moment coefficients are normalized by the reference
  **chord**, not the span. Scale by $c_{\text{ref}}/b_{\text{ref}}$ when
  comparing `Cl` or `Cn` against AVL output.

There is no run-case (`.run`) reader. Operating conditions are built in code as
a `FlightCondition`, and trim is set up through [Trim](./trim.md).

## Mass file

`MassReader` returns a `MassFileData`: unit multipliers, default gravity and
density, and a list of components.

Lines are classified in the same order as the legacy reader, which matters more
than it looks: a comment is stripped **before** classification, so an `=` or a
number inside a trailing comment cannot turn a real mass component into a unit
assignment.

| Line form | Meaning |
| --- | --- |
| `Lunit = 1.0` | length multiplier to SI; likewise `Munit`, `Tunit`, `g`, `rho` |
| `* sx sy sz` | scale subsequent component positions |
| `+ dx dy dz` | offset subsequent component positions |
| a leading letter | some other keyword line; skipped |
| seven or more numbers | a component |

A unit assignment with an unrecognized key falls through to component parsing
rather than consuming the line. A `*` or `+` directive needs three parseable
numbers and otherwise leaves the previous value untouched. A numeric run stops
at the first token that does not parse.

A component line is `mass x y z Ixx Iyy Izz [Ixy Ixz Iyz]`. The position has the
active scale and offset applied. The trailing comment becomes the component
name, which is how AVL mass files label their rows.

| Accessor | Returns |
| --- | --- |
| `total_mass()` | the summed mass, in SI |
| `center_of_gravity()` | the mass-weighted centroid, in SI |
| `total_properties()` | a `MassProperties`, in SI, with every component's inertia transferred to the combined CG by the parallel-axis theorem |

The three aggregate accessors **apply** the file's `Lunit` and `Munit`
multipliers and therefore return SI. The individual `MassComponent` records do
not: they stay in file units, with only the active `*` scale and `+` offset
applied to their positions. Mixing the two is the mistake to avoid.

If your geometry is not in SI, that conversion is yours to undo — the readers
never touch the geometry file's units.

Products of inertia are stored positive, matching
[Dynamic Modes](./dynamic-modes.md), and `MassProperties::inertia_tensor()`
applies the signs.

## Parse errors

A failed read returns `ZephyrError::Parse` carrying a `ParseLocation` — the
source path and the 1-based line number as separate fields, not baked into a
message string — plus the keyword being read. `read_path` attaches the path even
to errors raised by the text parser, so a message always names the file.

```rust
match AvlReader::new().read_path("aircraft.avl") {
    Ok(aircraft) => { let _ = aircraft; }
    Err(ZephyrError::Parse { location, message }) => {
        eprintln!("{location}: {message}");
    }
    Err(other) => eprintln!("{other}"),
}
```
