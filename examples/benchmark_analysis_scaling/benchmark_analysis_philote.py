# Single-analysis scaling benchmark: VLM via Philote gRPC.
#
# Sweeps spanwise panel count with fixed nx=6, runs mesh generation + solve
# multiple times per data point. Writes CSV to benchmark_scaling_out/.
#
# Prerequisites:
#   1. Start the mesh server:    ./mesh_philote_server localhost:50052
#   2. Start the VLM server:     ./vlm_philote_server localhost:50051

import sys
import time
import math
import numpy as np
import grpc
import openmdao.api as om
import philote_mdo.openmdao as pmom

span = 10.0
chord = 1.0
alpha_rad = 5.0 * math.pi / 180.0
nx = 6
n_repeats = 50
ny_values = [8, 12, 16, 24, 32, 48, 64, 96]

results = []

for ny in ny_values:
    n_stations = ny + 1
    total_panels = nx * 2 * ny

    prob = om.Problem()
    model = prob.model

    # Mesh discipline
    mesh_channel = grpc.insecure_channel("localhost:50052")
    model.add_subsystem("mesh_comp",
        pmom.RemoteExplicitComponent(channel=mesh_channel, nx=nx, ny=ny,
            cosine_chordwise_spacing="true", cosine_spanwise_spacing="false",
            y_duplicated="true"),
        promotes_inputs=["chord", "twist", "xle", "yle", "zle"],
        promotes_outputs=["reference_area", "reference_span",
                          "reference_chord"])

    # VLM discipline
    vlm_channel = grpc.insecure_channel("localhost:50051")
    model.add_subsystem("vlm",
        pmom.RemoteExplicitComponent(channel=vlm_channel,
            num_surfaces=1, surface_names="wing",
            surface_nx=str(nx), surface_ny=str(ny),
            surface_has_wake="true", surface_y_duplicated="true",
            surface_y_duplicate_plane="0.0", surface_component_index="0",
            surface_affected_by_freestream="true"),
        promotes_inputs=["alpha", "beta", "mach", "dynamic_pressure",
                         "roll_rate", "pitch_rate", "yaw_rate",
                         "reference_area", "reference_chord", "reference_span",
                         "moment_ref_x", "moment_ref_y", "moment_ref_z"],
        promotes_outputs=["CL", "CD", "CY", "Cl", "Cm", "Cn"])

    model.connect("mesh_comp.mesh", "vlm.wing")

    prob.setup()

    yle = np.linspace(0.0, span / 2.0, n_stations)
    prob["chord"] = np.full(n_stations, chord)
    prob["twist"] = np.zeros(n_stations)
    prob["xle"] = np.zeros(n_stations)
    prob["yle"] = yle
    prob["zle"] = np.zeros(n_stations)
    # reference_area, reference_chord, reference_span are now computed by
    # mesh_comp and connected automatically via promotes
    prob["moment_ref_x"] = 0.25 * chord
    prob["moment_ref_y"] = 0.0
    prob["moment_ref_z"] = 0.0
    prob["alpha"] = alpha_rad
    prob["beta"] = 0.0
    prob["mach"] = 0.0
    prob["dynamic_pressure"] = 0.5 * 1.225 * 25.0**2
    prob["roll_rate"] = 0.0
    prob["pitch_rate"] = 0.0
    prob["yaw_rate"] = 0.0

    # Warm-up
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

with open("benchmark_scaling_out/results_philote.csv", "w") as f:
    f.write("ny,panels,time_ms\n")
    for ny, panels, t in results:
        f.write(f"{ny},{panels},{t:.3f}\n")
print("Results written to benchmark_scaling_out/results_philote.csv")
