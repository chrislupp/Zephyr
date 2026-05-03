---
sidebar_position: 5
---

# MDO Integration

VLM can be used as a discipline in multidisciplinary design optimization (MDO) frameworks through the Philote-Cpp interface.

## Building with Philote Support

The Philote MDO interface requires additional dependencies:

- **Protobuf** — protocol buffer serialization
- **gRPC** — remote procedure call framework
- **PhiloteCpp** — Philote discipline interface library

Enable the Philote build option:

```bash
cmake -B build \
    -DCMAKE_BUILD_TYPE=Release \
    -DVLM_BUILD_PHILOTE=ON
```

This builds two server executables:

- `vlm_philote_server` — VLM aerodynamic analysis discipline
- `mesh_philote_server` — mesh generation discipline

## VLM Discipline Server

The VLM Philote discipline exposes aerodynamic analysis as a gRPC service. It accepts flight conditions and geometry parameters as inputs and returns aerodynamic coefficients as outputs.

### Inputs
- Flight conditions (angle of attack, sideslip, etc.)
- Geometry parameters (chord, twist, section positions)
- Mesh node positions (for aerostructural coupling)

### Outputs
- Aerodynamic force and moment coefficients
- Analytical adjoint gradients for all inputs

## Mesh Discipline Server

The mesh discipline provides parametric mesh generation as a standalone service, useful for coupling with structural solvers in aerostructural optimization.

## Usage with OpenMDAO

The Philote discipline interface is compatible with OpenMDAO and other MDO frameworks that support gRPC-based discipline evaluation. See the [Philote-Cpp documentation](https://github.com/chrislupp/Philote-Cpp) for framework integration details.
