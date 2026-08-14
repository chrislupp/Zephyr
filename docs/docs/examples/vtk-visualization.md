---
sidebar_position: 4
---

# VTK Visualization

The default `vtk` feature writes four dependency-free VTK XML products:
surface mesh, spanwise loads, wake, and velocity field.

```bash
cargo run --example vtk_visualization -- vtk_output
```

Both ASCII and binary files can be parsed back with `vtkio`. Exports use a
temporary sibling file and rename only after successful serialization.
