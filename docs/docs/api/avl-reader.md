---
sidebar_position: 15
---

# AVL and Mass Readers

`io::AvlReader` parses AVL geometry files into an `Aircraft`; `io::MassReader`
parses AVL mass files into a `MassFileData`. Both produce ordinary library
types — nothing downstream knows where the geometry came from.

The guide sheet, [AVL Input Files](../guide/avl-files.md), lists every supported
keyword and the deviations from AVL. This sheet is the type reference.

```rust
use zephyr_vlm::io::{AvlReader, MassReader};

let aircraft = AvlReader::new().read_path("aircraft.avl")?;
let mass = MassReader::read_path("aircraft.mass")?.total_properties();
```

## `AvlReader`

| Method | Reads |
| --- | --- |
| `read_path(path)` | a file, resolving relative `AFILE` paths against its directory |
| `read_str(content)` | text, with no base directory for airfoil files |
| `read_with_base(content, base)` | text, with an explicit base directory |
| `set_verbose(bool)` | print the configuration title and surface names to stderr |

Keywords are matched on their first four characters, case-insensitively.
Blank lines and lines starting with `#` or `!` are dropped before parsing.

A failed `AFILE` load is an **error**, not a silent fall back to a flat plate:
degrading the camber while reporting success would change the flow the solve
runs on. `read_path` attaches the source path to errors raised by the text
parser, so a message always names the file.

Surfaces, sections, airfoils, bodies, controls, and design-variable metadata are
all populated. `BODY` blocks are stored but **do not affect the solution**.

## `MassFileData`

| Field | Contents |
| --- | --- |
| `length_unit`, `mass_unit`, `time_unit` | multipliers to SI, as declared in the file |
| `gravity`, `density` | file defaults |
| `components` | the parsed `MassComponent`s |

| Accessor | Returns |
| --- | --- |
| `total_mass()` | the summed mass, in SI |
| `center_of_gravity()` | the mass-weighted centroid, in SI |
| `total_properties()` | a `MassProperties`, in SI, with each component's inertia transferred to the combined CG |

The three aggregate accessors **apply** the file's unit multipliers and return
SI. The individual `MassComponent` records do not: they stay in file units, with
only the active `*` scale and `+` offset applied to their positions. Mixing the
two is the mistake to avoid.

`MassComponent` carries `name` — taken from the row's trailing comment, which is
how AVL mass files label their components — `mass`, `position`, the three
principal inertias, and the three products of inertia, stored positive.

## `MassReader`

`read_path(path)` and `read_str(content)`. Both are associated functions; there
is no reader state.

[`AvlReader` rustdoc](https://docs.rs/zephyr-vlm/latest/zephyr_vlm/io/struct.AvlReader.html)
