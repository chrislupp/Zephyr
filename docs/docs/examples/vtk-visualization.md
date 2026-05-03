---
sidebar_position: 4
---

# VTK Visualization

This example demonstrates exporting VLM results to VTK format for visualization in ParaView.

## Requirements

VTK visualization requires building with VTK support:

```bash
cmake -B build -DVLM_ENABLE_VTK=ON
```

If VTK is not found at build time, the example still compiles but will report that VTK is unavailable at runtime.

## Overview

The program:
1. Creates a tapered wing with washout (twist variation)
2. Solves at a specified flight condition
3. Exports VTK files for surface mesh, spanwise loads, wake filaments, and velocity field

## Code

```cpp
#include "pegasus/vlm.h"

#include <iostream>

using namespace pegasus;

int main() {
    // Create a tapered wing with washout
    Surface wing("Wing");
    wing.SetChordwisePanelCount(8);
    wing.SetSpanwisePanelCount(24);
    wing.SetYDuplicate(0.0);

    Section root(Vec3(0.0, 0.0, 0.0), 2.0, 3.0);   // chord=2, twist=3 deg
    Section tip(Vec3(0.5, 6.0, 0.3), 0.8, -1.0);    // taper + washout
    wing.AddSection(root);
    wing.AddSection(tip);

    Aircraft aircraft("VTK Example");
    aircraft.AddSurface(wing);
    aircraft.SetReferenceArea(16.8);
    aircraft.SetReferenceChord(1.4);
    aircraft.SetReferenceSpan(12.0);

    VortexLatticeSolver solver;
    solver.SetAircraft(aircraft);

    FlightCondition condition;
    condition.alpha = 5.0 * kDegreesToRadians;
    SolutionState solution = solver.Solve(condition);

    // Export VTK files
    VtkExporter exporter;
    exporter.ExportSurfaceMesh(solver.GetMesh(), solution, "surface.vtp");
    exporter.ExportSpanwiseLoads(solver.GetMesh(), solution, "loads.vtp");
    exporter.ExportWake(solver.GetMesh(), "wake.vtp");

    std::cout << "VTK files exported." << std::endl;
    return 0;
}
```

## Running

```bash
./build/examples/vtk_visualization
```

## Viewing in ParaView

1. Open the exported `.vtp` files in ParaView
2. Apply a **Surface** representation for the mesh
3. Color by circulation ($\Gamma$) or pressure coefficient ($C_p$)
4. Use **Warp By Vector** filter to visualize deformed shapes
