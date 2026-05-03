# Copyright (c) 2024 VLM Authors
# SPDX-License-Identifier: AGPL-3.0-or-later
#
# OpenMDAO VLM analysis using Philote MDO disciplines.
#
# Prerequisites:
#   1. Start the mesh server:    ./mesh_philote_server localhost:50052
#   2. Start the VLM server:     ./vlm_philote_server localhost:50051
#   3. Run this script:          python openmdao_vlm_analysis.py

import numpy as np
import grpc
import openmdao.api as om
import philote_mdo.openmdao as pmom


# Wing geometry parameters
span = 10.0          # full span [m]
chord_val = 1.0      # constant chord [m]
alpha_deg = 5.0      # angle of attack [deg]
mach = 0.0           # incompressible
dynamic_pressure = 0.5 * 1.225 * 25.0**2  # ~383 Pa at 25 m/s

# Panel counts
nx = 6   # chordwise panels
ny = 12  # spanwise panels (half-span, will be y-duplicated)

n_stations = ny + 1

# Spanwise stations (half-span, linear spacing)
yle = np.linspace(0.0, span / 2.0, n_stations)

# Constant chord, zero twist, LE at x=0, z=0
chord_arr = np.full(n_stations, chord_val)
twist_arr = np.zeros(n_stations)
xle_arr = np.zeros(n_stations)
zle_arr = np.zeros(n_stations)

# Reference quantities
ref_area = span * chord_val      # planform area [m^2]
ref_chord = chord_val            # mean aerodynamic chord [m]
ref_span = span                  # full span [m]

# Build the OpenMDAO problem
prob = om.Problem()
model = prob.model

# Mesh generation discipline (port 50052)
mesh_channel = grpc.insecure_channel("localhost:50052")
model.add_subsystem(
    "mesh_comp",
    pmom.RemoteExplicitComponent(
        channel=mesh_channel,
        nx=nx,
        ny=ny,
        cosine_chordwise_spacing="true",
        cosine_spanwise_spacing="false",
    ),
    promotes_inputs=["chord", "twist", "xle", "yle", "zle"],
)

# VLM solver discipline (port 50051)
vlm_channel = grpc.insecure_channel("localhost:50051")
model.add_subsystem(
    "vlm",
    pmom.RemoteExplicitComponent(
        channel=vlm_channel,
        num_surfaces=1,
        surface_names="wing",
        surface_nx=str(nx),
        surface_ny=str(ny),
        surface_has_wake="true",
        surface_y_duplicated="true",
        surface_y_duplicate_plane="0.0",
        surface_component_index="0",
        surface_affected_by_freestream="true",
        reference_area=ref_area,
        reference_chord=ref_chord,
        reference_span=ref_span,
        moment_ref_x=0.25 * chord_val,
        moment_ref_y=0.0,
        moment_ref_z=0.0,
    ),
    promotes_inputs=["alpha", "beta", "mach", "dynamic_pressure",
                     "roll_rate", "pitch_rate", "yaw_rate"],
    promotes_outputs=["CL", "CD", "CY", "Cl", "Cm", "Cn"],
)

# Connect mesh output to VLM solver input
model.connect("mesh_comp.mesh", "vlm.wing")

prob.setup()

# Set mesh discipline inputs
prob["chord"] = chord_arr
prob["twist"] = twist_arr
prob["xle"] = xle_arr
prob["yle"] = yle
prob["zle"] = zle_arr

# Set flight condition inputs
prob["alpha"] = np.deg2rad(alpha_deg)
prob["beta"] = 0.0
prob["mach"] = mach
prob["dynamic_pressure"] = dynamic_pressure
prob["roll_rate"] = 0.0
prob["pitch_rate"] = 0.0
prob["yaw_rate"] = 0.0

# Run analysis
prob.run_model()

# Print results
print(f"\n{'='*50}")
print(f"VLM Analysis Results (alpha = {alpha_deg} deg)")
print(f"{'='*50}")
print(f"  CL = {prob['CL'][0]:.6f}")
print(f"  CD = {prob['CD'][0]:.6f}")
print(f"  CY = {prob['CY'][0]:.6f}")
print(f"  Cl = {prob['Cl'][0]:.6f}")
print(f"  Cm = {prob['Cm'][0]:.6f}")
print(f"  Cn = {prob['Cn'][0]:.6f}")
print(f"  L/D = {prob['CL'][0] / prob['CD'][0]:.2f}")
print(f"{'='*50}\n")
