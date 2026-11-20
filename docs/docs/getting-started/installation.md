---
sidebar_position: 1
---

# Installation

Zephyr requires Rust 1.93.1 or newer. Add the library to `Cargo.toml`:

```toml
[dependencies]
zephyr-vlm = { git = "https://github.com/chrislupp/Zephyr.git" }
```

VTK XML export is enabled by default and uses `vtkio`; no system VTK library
is required. For a minimal build:

```toml
zephyr-vlm = {
    git = "https://github.com/chrislupp/Zephyr.git",
    default-features = false,
}
```

Build and test a source checkout with:

```bash
cargo test --all-features
cargo test --no-default-features
```

The crate is Rust-only and forbids unsafe code. No CMake, C++ compiler,
Philote service, or Python environment is needed.
