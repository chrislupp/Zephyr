# Single-analysis scaling benchmark: OpenAeroStruct.
#
# Sweeps spanwise panel count with fixed nx=6, runs geometry + aero solve
# multiple times per data point for timing accuracy.

import time
import numpy as np
import openmdao.api as om
from openaerostruct.geometry.geometry_group import Geometry
from openaerostruct.aerodynamics.aero_groups import AeroPoint

span = 10.0
chord_init = 1.0
alpha_deg = 5.0
nx = 6
n_repeats = 50
ny_values = [8, 12, 16, 24, 32, 48, 64, 96]

import sys

results = []

for ny in ny_values:
    total_panels = nx * 2 * ny

    # Build mesh
    num_x = nx + 1
    num_y = ny + 1
    mesh = np.zeros((num_x, num_y, 3))
    y_stations = np.linspace(-span / 2.0, 0.0, num_y)
    for i in range(num_x):
        mesh[i, :, 0] = (i / nx) * chord_init
        mesh[i, :, 1] = y_stations

    surface = {
        "name": "wing",
        "symmetry": True,
        "S_ref_type": "projected",
        "mesh": mesh,
        "CL0": 0.0,
        "CD0": 0.0,
        "with_viscous": False,
        "with_wave": False,
        "k_lam": 0.05,
        "t_over_c_cp": np.array([0.12]),
        "c_max_t": 0.303,
    }

    # Warm-up run
    prob = om.Problem()
    model = prob.model
    indep = model.add_subsystem("indep", om.IndepVarComp(), promotes=["*"])
    indep.add_output("v", val=25.0, units="m/s")
    indep.add_output("alpha", val=alpha_deg, units="deg")
    indep.add_output("Mach_number", val=0.0)
    indep.add_output("re", val=1e6, units="1/m")
    indep.add_output("rho", val=1.225, units="kg/m**3")
    indep.add_output("cg", val=np.zeros(3), units="m")

    name = surface["name"]
    point_name = "aero_point_0"
    model.add_subsystem(name, Geometry(surface=surface))
    model.add_subsystem(point_name, AeroPoint(surfaces=[surface]),
                        promotes_inputs=["v", "alpha", "Mach_number", "re", "rho", "cg"])
    model.connect(name + ".mesh", point_name + "." + name + ".def_mesh")
    model.connect(name + ".mesh", point_name + ".aero_states." + name + "_def_mesh")
    model.connect(name + ".t_over_c", point_name + "." + name + "_perf.t_over_c")
    prob.setup()
    prob.run_model()

    # Timed runs
    t0 = time.perf_counter()
    for _ in range(n_repeats):
        prob.run_model()
    t_total = time.perf_counter() - t0
    per_run_ms = t_total / n_repeats * 1000

    results.append((ny, total_panels, per_run_ms))
    print(f"ny={ny} ({total_panels} panels): {per_run_ms:.2f} ms/run",
          file=sys.stderr)

# Write CSV to file
with open("benchmark_scaling_out/results_oas.csv", "w") as f:
    f.write("ny,panels,time_ms\n")
    for ny, panels, t in results:
        f.write(f"{ny},{panels},{t:.3f}\n")
print("Results written to benchmark_scaling_out/results_oas.csv")
