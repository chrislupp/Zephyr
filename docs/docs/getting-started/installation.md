---
sidebar_position: 1
---

# Installation

Zephyr is a library crate. There is no executable, no interactive shell, and no
input-file front end: an analysis is a Rust program that builds geometry, solves
it, and reads the result. Everything the library needs is a Rust toolchain.

## Toolchain

Rust 1.93.1 or newer. A source checkout pins it in `rust-toolchain.toml`, so a
`cargo` command run inside the repository selects the right compiler on its own.

The crate is `zephyr-vlm` on crates.io and is imported as `zephyr_vlm`.

```bash
cargo add zephyr-vlm
```

There is no C, C++, Fortran, or Python runtime dependency, no CMake step, and no
system VTK or BLAS to locate. The source forbids `unsafe`.

## Features

| Feature | Default | What it adds |
| --- | --- | --- |
| `vtk` | on | VTK XML export through the pure-Rust `vtkio` crate |
| `parallel` | off | Multithreaded influence assembly and geometry reverse sweeps, through `rayon` |
| `philote` | off | Philote-MDO gRPC disciplines, pulling in `tonic` and a `tokio` runtime |

A minimal build drops `vtkio`:

```bash
cargo add zephyr-vlm --no-default-features
```

With `vtk` off, the functions in `zephyr_vlm::visualization` return
`ZephyrError::FeatureDisabled` and `visualization::is_available()` returns
`false`. Nothing else changes.

`philote` is off by default on purpose. The gRPC stack and async runtime it
needs are not dependencies of the solver, and leaving the feature off keeps them
out of the dependency graph of every consumer that does not want MDO transport.
Nothing outside `zephyr_vlm::philote` and its example servers refers to it.

`parallel` is also off by default, so a standard build is single-threaded.
Results do not depend on the thread count: each source column writes a disjoint
slice of the influence operators, and the geometry reverse sweep splits rows
into a fixed number of chunks whose partial sums are accumulated in chunk order.
Thread count comes from `RAYON_NUM_THREADS`. Systems with fewer than 128
circulation unknowns stay on the serial loops, and the LU factorization is
serial either way.

```bash
cargo add zephyr-vlm --features parallel
```

## Following unreleased work

```toml
[dependencies]
zephyr-vlm = { git = "https://github.com/chrislupp/Zephyr.git" }
```

## Building a checkout

```bash
git clone https://github.com/chrislupp/Zephyr.git
cd Zephyr
cargo test --all-features
cargo test --no-default-features
```

The second command matters: it is what proves the library still builds and
passes without `vtkio`, `rayon`, and the Philote stack.

The checks the project runs before any change lands are:

```bash
cargo fmt --check
cargo clippy --all-targets --all-features -- -D warnings
cargo test --all-features
cargo test --release --all-features
cargo test --no-default-features
RUSTDOCFLAGS="-D warnings" cargo doc --all-features --no-deps
cargo package
```

Line coverage is enforced in CI with `cargo llvm-cov --fail-under-lines 92`.

## Running the examples

Every example in `examples/` is a complete program. The ones that need a feature
declare it, and Cargo refuses to build them without it.

```bash
cargo run --example basic_analysis
cargo run --release --example selective_analysis
cargo run --features philote --example philote_vlm_server
```

## Item-level reference

This site covers conventions, workflows, and the parts of the API where the
semantics matter more than the signature. Generated per-item documentation lives
on [docs.rs](https://docs.rs/zephyr-vlm), and `cargo doc --open` builds the same
pages for a checkout.
