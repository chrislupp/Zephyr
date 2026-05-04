# Copyright (c) 2024 VLM Authors
# SPDX-License-Identifier: AGPL-3.0-or-later
#
# OpenMDAO chord optimization using Philote MDO disciplines.
#
# Minimizes induced drag (CD) by varying the spanwise chord distribution
# subject to a lift constraint (CL = 0.5) and fixed planform area. The
# optimizer should converge toward an elliptic planform, which is the
# theoretical minimum-drag solution for a planar wing.
#
# Prerequisites:
#   1. Start the mesh server:    ./mesh_philote_server localhost:50052
#   2. Start the VLM server:     ./vlm_philote_server localhost:50051
#   3. Run this script:          python openmdao_vlm_optimization.py

import numpy as np
import grpc
import openmdao.api as om
import philote_mdo.openmdao as pmom


class AkimaChord(om.ExplicitComponent):
    """Interpolates chord from control points to spanwise stations using Akima.

    Control points and stations are defined in normalized span coordinates
    (eta in [0, 1]), so the interpolation is independent of halfspan.
    """

    def initialize(self):
        self.options.declare("n_cp", types=int)
        self.options.declare("n_stations", types=int)
        self.options.declare("eta_stations", types=np.ndarray)

    def setup(self):
        n_cp = self.options["n_cp"]
        n_st = self.options["n_stations"]
        self.add_input("chord_cp", shape=(n_cp,), units="m")
        self.add_output("chord", shape=(n_st,), units="m")

    def setup_partials(self):
        self.declare_partials("chord", "chord_cp")

    def compute(self, inputs, outputs):
        from scipy.interpolate import Akima1DInterpolator
        n_cp = self.options["n_cp"]
        eta_st = self.options["eta_stations"]

        eta_cp = np.linspace(0.0, 1.0, n_cp)
        spline = Akima1DInterpolator(eta_cp, inputs["chord_cp"])
        outputs["chord"] = spline(eta_st)

    def compute_partials(self, inputs, partials):
        from scipy.interpolate import Akima1DInterpolator
        n_cp = self.options["n_cp"]
        n_st = self.options["n_stations"]
        eta_st = self.options["eta_stations"]

        eta_cp = np.linspace(0.0, 1.0, n_cp)
        cp = inputs["chord_cp"].ravel()

        # Compute Jacobian via finite differences on the spline
        jac = np.zeros((n_st, n_cp))
        h = 1e-8
        for k in range(n_cp):
            cp_p = cp.copy()
            cp_p[k] += h
            sp_p = Akima1DInterpolator(eta_cp, cp_p)
            cp_m = cp.copy()
            cp_m[k] -= h
            sp_m = Akima1DInterpolator(eta_cp, cp_m)
            jac[:, k] = (sp_p(eta_st) - sp_m(eta_st)) / (2.0 * h)

        partials["chord", "chord_cp"] = jac


# Wing geometry parameters
halfspan_init = 5.0  # initial half-span [m]
chord_init = 1.0     # initial chord [m]
alpha_deg = 5.0      # angle of attack [deg]
mach = 0.0           # incompressible
dynamic_pressure = 0.5 * 1.225 * 25.0**2  # ~383 Pa at 25 m/s

# Target lift coefficient
CL_target = 0.5

# Panel counts
nx = 6   # chordwise panels
ny = 12  # spanwise panels (half-span, y-duplicated)

n_stations = ny + 1
n_cp = 5  # chord B-spline control points

# Normalized spanwise stations (eta in [0, 1])
eta_stations = np.linspace(0.0, 1.0, n_stations)

# Initial geometry: constant chord, zero twist
chord_cp_init = np.full(n_cp, chord_init)
twist_arr = np.zeros(n_stations)
xle_arr = np.zeros(n_stations)
zle_arr = np.zeros(n_stations)

# Initial planform area (used as constraint target)
ref_area = 2.0 * halfspan_init * chord_init

# Build the OpenMDAO problem
prob = om.Problem()
model = prob.model

# Chord interpolation (control points -> stations)
model.add_subsystem(
    "chord_interp",
    AkimaChord(n_cp=n_cp, n_stations=n_stations, eta_stations=eta_stations),
    promotes_inputs=[("chord_cp", "chord_cp")],
    promotes_outputs=["chord"],
)

# Mesh generation discipline (port 50052)
mesh_channel = grpc.insecure_channel("localhost:50052")
model.add_subsystem(
    "mesh_comp",
    pmom.RemoteExplicitComponent(
        channel=mesh_channel,
        nx=nx,
        ny=ny,
        cosine_chordwise_spacing="true",
        y_duplicated="true",
    ),
    promotes_inputs=["chord", "twist", "xle", "zle", "halfspan", "eta"],
    promotes_outputs=["reference_area", "reference_span", "reference_chord"],
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
    ),
    promotes_inputs=["alpha", "beta", "mach", "dynamic_pressure",
                     "roll_rate", "pitch_rate", "yaw_rate",
                     "reference_area", "reference_chord", "reference_span",
                     "moment_ref_x", "moment_ref_y", "moment_ref_z"],
    promotes_outputs=["CL", "CD", "CY", "Cl", "Cm", "Cn"],
)

# Connect mesh output to VLM solver input
model.connect("mesh_comp.mesh", "vlm.wing")

# Optimizer
prob.driver = om.ScipyOptimizeDriver()
prob.driver.options["optimizer"] = "SLSQP"
prob.driver.options["tol"] = 1e-8
prob.driver.options["maxiter"] = 200

# Record optimization iterations to sqlite
recorder = om.SqliteRecorder("openmdao_vlm_optimization_out/opt_history.sql")
prob.driver.add_recorder(recorder)
prob.driver.recording_options["includes"] = ["CD"]
prob.driver.recording_options["record_derivatives"] = True

# Design variable: chord control points
model.add_design_var("chord_cp", lower=0.1, upper=3.0)

# Design variable: half-span (upper bound: 50% growth from initial)
model.add_design_var("halfspan", lower=2.0, upper=7.5)

# Design variable: angle of attack (free to meet CL constraint)
model.add_design_var("alpha", lower=np.deg2rad(-10.0), upper=np.deg2rad(15.0))

# Objective: minimize drag
model.add_objective("CD")

# Constraints
model.add_constraint("CL", equals=CL_target)
model.add_constraint("reference_area", equals=ref_area)

prob.setup()

# Set initial values
prob["chord_cp"] = chord_cp_init
prob["twist"] = twist_arr
prob["xle"] = xle_arr
prob["zle"] = zle_arr
prob["halfspan"] = halfspan_init
prob["eta"] = eta_stations

# reference_area, reference_chord, reference_span computed by mesh_comp
prob["moment_ref_x"] = 0.25 * chord_init
prob["moment_ref_y"] = 0.0
prob["moment_ref_z"] = 0.0

# Set flight condition inputs
prob["alpha"] = np.deg2rad(alpha_deg)
prob["beta"] = 0.0
prob["mach"] = mach
prob["dynamic_pressure"] = dynamic_pressure
prob["roll_rate"] = 0.0
prob["pitch_rate"] = 0.0
prob["yaw_rate"] = 0.0

# Run optimization
prob.run_driver()
prob.cleanup()

# Print results
opt_chord = prob["chord"]
opt_alpha = np.rad2deg(prob["alpha"][0])
opt_halfspan = prob["halfspan"][0]
opt_span = 2.0 * opt_halfspan
opt_yle = eta_stations * opt_halfspan

print(f"\n{'='*50}")
print(f"Optimization Results")
print(f"{'='*50}")
print(f"  CL = {prob['CL'][0]:.6f}  (target: {CL_target})")
print(f"  CD = {prob['CD'][0]:.6f}")
print(f"  L/D = {prob['CL'][0] / prob['CD'][0]:.2f}")
print(f"  alpha = {opt_alpha:.3f} deg")
print(f"  half-span = {opt_halfspan:.4f} m  (initial: {halfspan_init})")
print(f"  full span = {opt_span:.4f} m")
print(f"  planform area = {prob['reference_area'][0]:.4f} m^2  (target: {ref_area})")
print(f"\n  Chord control points: {prob['chord_cp'].ravel()}")
print(f"\n  Chord distribution (root to tip):")
for i, (y, c) in enumerate(zip(opt_yle, opt_chord)):
    print(f"    y={y:6.3f} m  chord={c:.4f} m")
print(f"{'='*50}\n")

# Compare to elliptic planform
c_root = opt_chord[0]
elliptic = c_root * np.sqrt(np.maximum(1.0 - (2.0 * opt_yle / opt_span) ** 2, 0.0))
print("Comparison to elliptic planform (c_root matched):")
print(f"  {'y [m]':>8s}  {'optimized':>10s}  {'elliptic':>10s}  {'ratio':>8s}")
for y, c_opt, c_ell in zip(opt_yle, opt_chord, elliptic):
    ratio = c_opt / c_ell if c_ell > 1e-6 else float("inf")
    print(f"  {y:8.3f}  {c_opt:10.4f}  {c_ell:10.4f}  {ratio:8.3f}")
print()

# ── Plotting ──────────────────────────────────────────────────────────────
import matplotlib.pyplot as plt

# Full-span arrays (mirror the half-span data)
y_full = np.concatenate([-opt_yle[::-1], opt_yle[1:]])
chord_full = np.concatenate([opt_chord[::-1], opt_chord[1:]])
elliptic_full = np.concatenate([elliptic[::-1], elliptic[1:]])

subtitle = (
    f"$C_L$={prob['CL'][0]:.4f}, $C_D$={prob['CD'][0]:.6f}, "
    f"$L/D$={prob['CL'][0]/prob['CD'][0]:.1f}, "
    f"$\\alpha$={opt_alpha:.2f}$^\\circ$, "
    f"span={opt_span:.2f} m"
)

# --- Figure 1: Chord distribution ---
fig1, ax1 = plt.subplots(figsize=(7, 5))
y_ellipse_fine = np.linspace(-opt_halfspan, opt_halfspan, 200)
ellipse_fine = c_root * np.sqrt(np.maximum(1.0 - (y_ellipse_fine / opt_halfspan) ** 2, 0.0))
ax1.plot(y_ellipse_fine, ellipse_fine, "k--", linewidth=1.0, label="Elliptic (matched $c_\\mathrm{root}$)")
ax1.plot(y_full, chord_full, "o-", color="tab:blue", markersize=4, linewidth=1.5, label="Optimized")
ax1.set_xlabel("Spanwise station $y$ [m]")
ax1.set_ylabel("Chord [m]")
ax1.set_title(f"Chord Distribution\n{subtitle}")
ax1.legend(loc="lower center")
ax1.set_xlim(-opt_halfspan * 1.05, opt_halfspan * 1.05)
ax1.set_ylim(bottom=0)
ax1.grid(True, alpha=0.3)
fig1.tight_layout()
fig1.savefig("openmdao_vlm_optimization_out/chord_distribution.png", dpi=150)

# --- Figure 2: Wing planform ---
fig2, ax2 = plt.subplots(figsize=(10, 4))
xle_full = np.concatenate([xle_arr[::-1], xle_arr[1:]])
le_x = xle_full - 0.25 * chord_full
te_x = xle_full + 0.75 * chord_full

# Optimized planform (y horizontal, x vertical with flow going left-to-right)
ax2.fill_between(y_full, -le_x, -te_x, alpha=0.25, color="tab:blue", label="Optimized")
ax2.plot(y_full, -le_x, "-", color="tab:blue", linewidth=1.5)
ax2.plot(y_full, -te_x, "-", color="tab:blue", linewidth=1.5)
ax2.plot([y_full[0], y_full[0]], [-le_x[0], -te_x[0]], "-", color="tab:blue", linewidth=1.5)
ax2.plot([y_full[-1], y_full[-1]], [-le_x[-1], -te_x[-1]], "-", color="tab:blue", linewidth=1.5)

# Elliptic planform outline (smooth)
y_ell_fine = np.linspace(-opt_halfspan, opt_halfspan, 200)
c_ell_fine = c_root * np.sqrt(np.maximum(1.0 - (y_ell_fine / opt_halfspan) ** 2, 0.0))
le_ell_fine = -0.25 * c_ell_fine
te_ell_fine = 0.75 * c_ell_fine
ax2.plot(y_ell_fine, -le_ell_fine, "k--", linewidth=1.0, label="Elliptic")
ax2.plot(y_ell_fine, -te_ell_fine, "k--", linewidth=1.0)

ax2.set_xlabel("$y$ [m]")
ax2.set_ylabel("$x$ [m]")
ax2.set_title(f"Wing Planform\n{subtitle}")
ax2.set_aspect("equal")
ax2.grid(True, alpha=0.3)
fig2.tight_layout()
fig2.savefig("openmdao_vlm_optimization_out/wing_planform.png", dpi=150)

# --- Figure 3: Optimization convergence history ---
import sqlite3

db_path = "openmdao_vlm_optimization_out/opt_history.sql"
conn = sqlite3.connect(db_path)

# Get iteration coordinates that have derivative evaluations (major iterations)
deriv_coords = {
    row[0]
    for row in conn.execute(
        "SELECT iteration_coordinate FROM driver_derivatives"
    ).fetchall()
}

# Get all driver iterations and filter to major ones
rows = conn.execute(
    "SELECT iteration_coordinate, outputs FROM driver_iterations"
    " ORDER BY counter"
).fetchall()
conn.close()

import json

major_cd = []
for coord, outputs_blob in rows:
    if coord in deriv_coords:
        outputs = json.loads(outputs_blob)
        # Find the CD key (may be promoted or absolute path)
        cd_val = None
        for key, val in outputs.items():
            if key.endswith(".CD") or key == "CD":
                cd_val = val[0] if isinstance(val, list) else val
                break
        if cd_val is not None:
            major_cd.append(cd_val)

fig3, ax3 = plt.subplots(figsize=(7, 5))
ax3.plot(range(1, len(major_cd) + 1), major_cd, "o-", color="tab:red",
         markersize=4, linewidth=1.5)
ax3.set_xlabel("Major Iteration")
ax3.set_ylabel("$C_D$")
ax3.set_title("Optimization Convergence History")
ax3.grid(True, alpha=0.3)
fig3.tight_layout()
fig3.savefig("openmdao_vlm_optimization_out/convergence_history.png", dpi=150)

plt.show()
print("Plots saved to openmdao_vlm_optimization_out/")
