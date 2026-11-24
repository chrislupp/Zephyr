---
sidebar_position: 17
---

# Errors

Fallible calls return `zephyr_vlm::Result<T>`, an alias for
`Result<T, ZephyrError>`. `ZephyrError` implements `std::error::Error` through
`thiserror`, so it composes with `anyhow`, `?`, and the rest of the ecosystem
normally.

```rust
use zephyr_vlm::{Result, ZephyrError, io::AvlReader};

fn load(path: &str) -> Result<()> {
    match AvlReader::new().read_path(path) {
        Ok(aircraft) => { let _ = aircraft; Ok(()) }
        Err(ZephyrError::Parse { location, message }) => {
            eprintln!("{location}: {message}");
            Ok(())
        }
        Err(other) => Err(other),
    }
}
```

## Variants

| Variant | Raised when |
| --- | --- |
| `Io { path, source }` | a file could not be read or written; the path is kept |
| `Parse { location, message }` | text input was invalid; `location` carries the path and line |
| `InvalidGeometry(String)` | an aircraft, airfoil, or adjoint configuration is unusable |
| `InvalidMesh(String)` | a direct mesh is invalid, or does not describe the solved lattice |
| `InvalidFlightCondition(String)` | a non-finite or out-of-range condition, or an unknown control name |
| `InvalidTrim(String)` | a trim problem is misconfigured — not square, or a control not owned |
| `TrimFailed(String)` | a trim solve did not converge |
| `SingularSystem` | the aerodynamic system is singular or ill-conditioned |
| `Eigensolver(String)` | the eigensolver failed, or its inputs were rejected |
| `InvalidSolution` | a valid, matching aerodynamic solution was required |
| `FeatureDisabled(&'static str)` | a call needs a Cargo feature this build does not have |

## `ParseLocation`

```rust
// `location` is the ParseLocation from a ZephyrError::Parse.
if let (Some(path), Some(line)) = (&location.path, location.line) {
    eprintln!("{}:{line}", path.display());
}
```

The path and the 1-based line number are **separate fields**, not baked into a
rendered message, so a caller can jump to the line or group errors by file
without parsing a string back apart. `Display` renders them for a message.

Constructors, for code building on the readers: `ZephyrError::io`,
`ZephyrError::parse`, `ZephyrError::parse_path`, and `with_source_path`, which
attaches a path to an error raised by a text parser.

## Things that are not errors

Two outcomes are deliberately **not** errors, and both are read wrongly often
enough to be worth stating:

- **`TrimSolver` nonconvergence** is a `TrimStatus` on a returned `TrimResult`,
  which carries its last state either way. Check `status` before the numbers.
  (`MeshTrimSolver` does the opposite: a failure there *is* an error.)
- **An invalid force result** is `AerodynamicResults` with
  `is_valid == false` and every field at its default. Check `is_valid`; a
  default result is all zeros, which is a plausible-looking lie.

## Mach and validity

`VortexLatticeSolver::solve` rejects Mach at or above 1 with
`InvalidFlightCondition` rather than returning the incompressible answer. It
also rejects non-finite angles, rates, Mach, and dynamic pressure, and negative
Mach or dynamic pressure. A NaN that reached the right-hand side would otherwise
produce an all-NaN circulation reported as a valid solution.

[`ZephyrError` rustdoc](https://docs.rs/zephyr-vlm/latest/zephyr_vlm/error/enum.ZephyrError.html)
