# VLM - Vortex Lattice Method Library

[![Build and Tests](https://github.com/chrislupp/Pegasus/actions/workflows/build.yml/badge.svg?branch=main)](https://github.com/chrislupp/Pegasus/actions/workflows/build.yml)
[![codecov](https://codecov.io/gh/chrislupp/Pegasus/branch/main/graph/badge.svg)](https://codecov.io/gh/chrislupp/Pegasus)

A modern C++17 implementation of the Vortex Lattice Method for aerodynamic analysis, ported from OptVL with enhancements for automatic differentiation and MDO integration.

## Features

- **Core VLM Solver**: Horseshoe vortex elements with Biot-Savart influence computation
- **Force Computation**: Near-field Kutta-Joukowski and far-field Trefftz plane drag
- **Stability Analysis**: Complete stability and body-axis derivatives
- **Trim Solver**: Newton-Raphson iteration for constrained flight conditions
- **Eigenmode Analysis**: Longitudinal and lateral-directional dynamic modes
- **Automatic Differentiation**: Forward-mode (dual numbers) and reverse-mode (tape-based)
- **AVL Compatibility**: Parse AVL geometry and mass property files

## Requirements

- CMake 3.16+
- C++17 compiler (GCC 8+, Clang 7+, MSVC 2019+)
- Eigen 5.0+ (automatically found or specify `CMAKE_PREFIX_PATH`)

## Building

```bash
mkdir build && cd build
cmake ..
make -j4
```

### Build Options

| Option | Default | Description |
|--------|---------|-------------|
| `VLM_BUILD_TESTS` | ON | Build unit and integration tests |
| `VLM_BUILD_EXAMPLES` | ON | Build example programs |
| `VLM_ENABLE_OPENMP` | ON | Enable OpenMP parallelization |

```bash
# Example: Release build without tests
cmake -DCMAKE_BUILD_TYPE=Release -DVLM_BUILD_TESTS=OFF ..
```

### Running Tests

```bash
cd build
ctest --output-on-failure
# Or run directly:
./tests/vlm_unit_tests
./tests/vlm_integration_tests
```

## Quick Start

### Basic Analysis

```cpp
#include "vlm/vlm.h"

using namespace vlm;

int main() {
    // Create a rectangular wing
    Surface wing("Wing");
    wing.SetChordwisePanelCount(6);
    wing.SetSpanwisePanelCount(20);
    wing.SetYDuplicate(0.0);  // Mirror about Y=0

    Section root(Vec3(0.0, 0.0, 0.0), 1.0);  // chord = 1.0
    Section tip(Vec3(0.0, 5.0, 0.0), 1.0);   // semi-span = 5.0
    wing.AddSection(root);
    wing.AddSection(tip);

    // Create aircraft
    Aircraft aircraft("Test Wing");
    aircraft.AddSurface(wing);
    aircraft.SetReferenceArea(10.0);
    aircraft.SetReferenceChord(1.0);
    aircraft.SetReferenceSpan(10.0);

    // Solve
    VortexLatticeSolver solver;
    solver.SetAircraft(aircraft);

    FlightCondition condition;
    condition.alpha = 5.0 * kDegreesToRadians;

    SolutionState solution = solver.Solve(condition);

    // Compute forces
    ForceCalculator calculator;
    AerodynamicResults results = calculator.ComputeForces(
        aircraft, solver.GetMesh(), solution, condition);

    std::cout << "CL = " << results.CL << std::endl;
    std::cout << "CD = " << results.CD << std::endl;
    std::cout << "L/D = " << results.CL / results.CD << std::endl;

    return 0;
}
```

### Loading AVL Files

```cpp
#include "vlm/vlm.h"

using namespace vlm;

int main() {
    AvlReader reader;
    Aircraft aircraft = reader.ReadFile("aircraft.avl");

    VortexLatticeSolver solver;
    solver.SetAircraft(aircraft);

    FlightCondition condition;
    condition.alpha = 5.0 * kDegreesToRadians;

    SolutionState solution = solver.Solve(condition);

    ForceCalculator calculator;
    AerodynamicResults results = calculator.ComputeForces(
        aircraft, solver.GetMesh(), solution, condition);

    std::cout << "Aircraft: " << aircraft.GetTitle() << std::endl;
    std::cout << "CL = " << results.CL << std::endl;

    return 0;
}
```

### Stability Derivatives

```cpp
#include "vlm/vlm.h"

using namespace vlm;

int main() {
    // ... set up aircraft and solver ...

    FlightCondition condition;
    condition.alpha = 5.0 * kDegreesToRadians;
    solver.Solve(condition);

    StabilityDerivativeCalculator deriv_calc;
    StabilityDerivatives derivs = deriv_calc.ComputeDerivatives(
        aircraft, solver, condition);

    std::cout << "CL_alpha = " << derivs.CL_alpha << " /rad" << std::endl;
    std::cout << "Cm_alpha = " << derivs.Cm_alpha << " /rad" << std::endl;
    std::cout << "Cl_p = " << derivs.Cl_p << " (roll damping)" << std::endl;

    return 0;
}
```

### Trim Analysis

```cpp
#include "vlm/vlm.h"

using namespace vlm;

int main() {
    // ... set up aircraft and solver ...

    TrimSolver trim_solver;
    TrimConstraints constraints = TrimPresets::LevelFlight(0.5);  // Target CL = 0.5

    TrimResult result = trim_solver.Solve(aircraft, solver, constraints);

    if (result.converged) {
        std::cout << "Trimmed alpha = "
                  << result.state.alpha * kRadiansToDegrees << " deg" << std::endl;
        std::cout << "CL = " << result.aero.CL << std::endl;
    }

    return 0;
}
```

### Eigenmode Analysis

```cpp
#include "vlm/vlm.h"

using namespace vlm;

int main() {
    // ... compute stability derivatives ...

    MassProperties mass;
    mass.mass = 1000.0;  // kg
    mass.Ixx = 1000.0;   // kg*m^2
    mass.Iyy = 2000.0;
    mass.Izz = 2500.0;

    FlightState flight;
    flight.velocity = 50.0;  // m/s
    flight.density = 1.225;  // kg/m^3

    EigenmodeAnalysis analysis;
    EigenmodeResults results = analysis.AnalyzeLongitudinal(
        aircraft, mass, derivs, flight);

    for (const auto& mode : results.modes) {
        std::cout << GetModeTypeName(mode.type) << ": ";
        if (mode.IsOscillatory()) {
            std::cout << "freq=" << mode.GetFrequency() << " rad/s, "
                      << "damping=" << mode.GetDampingRatio();
        }
        std::cout << (mode.IsStable() ? " [STABLE]" : " [UNSTABLE]") << std::endl;
    }

    return 0;
}
```

## Project Structure

```
VLM/
├── include/vlm/
│   ├── vlm.h                 # Main facade header
│   ├── core/                 # Types, constants
│   ├── geometry/             # Airfoil, Section, Surface, Aircraft
│   ├── io/                   # AVL/mass file readers
│   ├── solver/               # VLM solver, influence matrix
│   ├── aero/                 # Forces, stability derivatives
│   ├── trim/                 # Trim solver
│   ├── dynamics/             # Eigenmode analysis
│   ├── ad/                   # Automatic differentiation
│   └── utils/                # Cubic spline, math utilities
├── src/                      # Implementation files
├── tests/                    # Unit and integration tests
├── examples/                 # Example programs
└── CMakeLists.txt
```

## Examples

Three example programs are provided:

- **basic_analysis**: Angle of attack sweep for a rectangular wing
- **trim_flight**: Trim solver for level flight at various CL targets
- **stability_analysis**: Full stability derivatives and eigenmode analysis

Run examples from the build directory:

```bash
./examples/basic_analysis
./examples/trim_flight
./examples/stability_analysis
```

## API Reference

### Core Classes

| Class | Description |
|-------|-------------|
| `Aircraft` | Complete aircraft configuration |
| `Surface` | Lifting surface definition |
| `Section` | Wing section (chord, twist, airfoil) |
| `VortexLatticeSolver` | Main VLM solver |
| `ForceCalculator` | Force and moment computation |
| `StabilityDerivativeCalculator` | Stability derivatives |
| `TrimSolver` | Newton-Raphson trim |
| `EigenmodeAnalysis` | Dynamic stability modes |
| `AvlReader` | AVL file parser |
| `MassReader` | Mass property file parser |

### Key Structures

| Structure | Description |
|-----------|-------------|
| `FlightCondition` | Alpha, beta, rates, Mach |
| `SolutionState` | Circulation distribution |
| `AerodynamicResults` | Forces, moments, coefficients |
| `StabilityDerivatives` | All stability derivatives |
| `TrimResult` | Trim solution state |
| `EigenmodeResults` | Dynamic mode eigenvalues |

## Documentation

Full documentation is available at [https://chrislupp.github.io/Pegasus/](https://chrislupp.github.io/Pegasus/).

## License

GPL-3.0-or-later

## Acknowledgments

Ported from OptVL (Python/Fortran Vortex Lattice Method implementation).
