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


class PlanformArea(om.ExplicitComponent):
    """Computes planform area from chord and yle via trapezoidal integration.

    For a y-duplicated (half-span) wing, multiplies by 2 for full-span area.
    """

    def initialize(self):
        self.options.declare("n_stations", types=int)

    def setup(self):
        n = self.options["n_stations"]
        self.add_input("chord", shape=(n,), units="m")
        self.add_input("yle", shape=(n,), units="m")
        self.add_output("planform_area", shape=(1,), units="m**2")

    def setup_partials(self):
        n = self.options["n_stations"]
        self.declare_partials("planform_area", "chord",
                              rows=np.zeros(n, dtype=int),
                              cols=np.arange(n))
        self.declare_partials("planform_area", "yle",
                              rows=np.zeros(n, dtype=int),
                              cols=np.arange(n))

    def compute(self, inputs, outputs):
        c = inputs["chord"]
        y = inputs["yle"]
        dy = np.diff(y)
        outputs["planform_area"] = 2.0 * np.sum(0.5 * (c[:-1] + c[1:]) * dy)

    def compute_partials(self, inputs, partials):
        c = inputs["chord"]
        y = inputs["yle"]
        dy = np.diff(y)
        n = len(c)

        da_dc = np.zeros(n)
        for j in range(n):
            if j > 0:
                da_dc[j] += dy[j - 1]
            if j < n - 1:
                da_dc[j] += dy[j]
        partials["planform_area", "chord"] = da_dc

        da_dy = np.zeros(n)
        for j in range(n - 1):
            avg_c = 0.5 * (c[j] + c[j + 1])
            da_dy[j] -= 2.0 * avg_c
            da_dy[j + 1] += 2.0 * avg_c
        partials["planform_area", "yle"] = da_dy


class AkimaChord(om.ExplicitComponent):
    """Interpolates chord from control points to spanwise stations using Akima.

    Uses scipy's Akima1DInterpolator.
    """

    def initialize(self):
        self.options.declare("n_cp", types=int)
        self.options.declare("n_stations", types=int)
        self.options.declare("y_stations", types=np.ndarray)

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
        y_st = self.options["y_stations"]
        half_span = y_st[-1]

        y_cp = np.linspace(0.0, half_span, n_cp)
        spline = Akima1DInterpolator(y_cp, inputs["chord_cp"])
        outputs["chord"] = spline(y_st)

    def compute_partials(self, inputs, partials):
        from scipy.interpolate import Akima1DInterpolator
        n_cp = self.options["n_cp"]
        n_st = self.options["n_stations"]
        y_st = self.options["y_stations"]
        half_span = y_st[-1]

        y_cp = np.linspace(0.0, half_span, n_cp)
        cp = inputs["chord_cp"].ravel()

        # Compute Jacobian via finite differences on the spline
        jac = np.zeros((n_st, n_cp))
        h = 1e-8
        for k in range(n_cp):
            cp_p = cp.copy()
            cp_p[k] += h
            sp_p = Akima1DInterpolator(y_cp, cp_p)
            cp_m = cp.copy()
            cp_m[k] -= h
            sp_m = Akima1DInterpolator(y_cp, cp_m)
            jac[:, k] = (sp_p(y_st) - sp_m(y_st)) / (2.0 * h)

        partials["chord", "chord_cp"] = jac


# Wing geometry parameters
span = 10.0          # full span [m]
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

# Spanwise stations (half-span, linear spacing)
yle = np.linspace(0.0, span / 2.0, n_stations)

# Initial geometry: constant chord, zero twist
chord_cp_init = np.full(n_cp, chord_init)
twist_arr = np.zeros(n_stations)
xle_arr = np.zeros(n_stations)
zle_arr = np.zeros(n_stations)

# Reference quantities
ref_area = span * chord_init
ref_chord = chord_init
ref_span = span

# Build the OpenMDAO problem
prob = om.Problem()
model = prob.model

# Chord interpolation (control points -> stations)
model.add_subsystem(
    "chord_interp",
    AkimaChord(n_cp=n_cp, n_stations=n_stations, y_stations=yle),
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
    ),
    promotes_inputs=["alpha", "beta", "mach", "dynamic_pressure",
                     "roll_rate", "pitch_rate", "yaw_rate",
                     "reference_area", "reference_chord", "reference_span",
                     "moment_ref_x", "moment_ref_y", "moment_ref_z"],
    promotes_outputs=["CL", "CD", "CY", "Cl", "Cm", "Cn"],
)

# Planform area computation
model.add_subsystem(
    "area_comp",
    PlanformArea(n_stations=n_stations),
    promotes_inputs=["chord", "yle"],
    promotes_outputs=["planform_area"],
)

# Connect mesh output to VLM solver input
model.connect("mesh_comp.mesh", "vlm.wing")

# Optimizer
prob.driver = om.ScipyOptimizeDriver()
prob.driver.options["optimizer"] = "SLSQP"
prob.driver.options["tol"] = 1e-8
prob.driver.options["maxiter"] = 200

# Design variable: chord control points
model.add_design_var("chord_cp", lower=0.1, upper=3.0)

# Design variable: angle of attack (free to meet CL constraint)
model.add_design_var("alpha", lower=np.deg2rad(-10.0), upper=np.deg2rad(15.0))

# Objective: minimize drag
model.add_objective("CD")

# Constraints
model.add_constraint("CL", equals=CL_target)
model.add_constraint("planform_area", equals=ref_area)

prob.setup()

# Set initial values
prob["chord_cp"] = chord_cp_init
prob["twist"] = twist_arr
prob["xle"] = xle_arr
prob["yle"] = yle
prob["zle"] = zle_arr

# Set reference quantities
prob["reference_area"] = ref_area
prob["reference_chord"] = ref_chord
prob["reference_span"] = ref_span
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

# Print results
opt_chord = prob["chord"]
opt_alpha = np.rad2deg(prob["alpha"][0])

print(f"\n{'='*50}")
print(f"Optimization Results")
print(f"{'='*50}")
print(f"  CL = {prob['CL'][0]:.6f}  (target: {CL_target})")
print(f"  CD = {prob['CD'][0]:.6f}")
print(f"  L/D = {prob['CL'][0] / prob['CD'][0]:.2f}")
print(f"  alpha = {opt_alpha:.3f} deg")
print(f"  planform area = {prob['planform_area'][0]:.4f} m^2  (target: {ref_area})")
print(f"\n  Chord control points: {prob['chord_cp'].ravel()}")
print(f"\n  Chord distribution (root to tip):")
for i, (y, c) in enumerate(zip(yle, opt_chord)):
    print(f"    y={y:6.3f} m  chord={c:.4f} m")
print(f"{'='*50}\n")

# Compare to elliptic planform
c_root = opt_chord[0]
elliptic = c_root * np.sqrt(np.maximum(1.0 - (2.0 * yle / span) ** 2, 0.0))
print("Comparison to elliptic planform (c_root matched):")
print(f"  {'y [m]':>8s}  {'optimized':>10s}  {'elliptic':>10s}  {'ratio':>8s}")
for y, c_opt, c_ell in zip(yle, opt_chord, elliptic):
    ratio = c_opt / c_ell if c_ell > 1e-6 else float("inf")
    print(f"  {y:8.3f}  {c_opt:10.4f}  {c_ell:10.4f}  {ratio:8.3f}")
print()

# ── Plotting ──────────────────────────────────────────────────────────────
import matplotlib.pyplot as plt

# Full-span arrays (mirror the half-span data)
y_full = np.concatenate([-yle[::-1], yle[1:]])
chord_full = np.concatenate([opt_chord[::-1], opt_chord[1:]])
elliptic_full = np.concatenate([elliptic[::-1], elliptic[1:]])

subtitle = (
    f"$C_L$={prob['CL'][0]:.4f}, $C_D$={prob['CD'][0]:.6f}, "
    f"$L/D$={prob['CL'][0]/prob['CD'][0]:.1f}, "
    f"$\\alpha$={opt_alpha:.2f}$^\\circ$"
)

# --- Figure 1: Chord distribution ---
fig1, ax1 = plt.subplots(figsize=(7, 5))
y_ellipse_fine = np.linspace(-span / 2, span / 2, 200)
ellipse_fine = c_root * np.sqrt(np.maximum(1.0 - (2.0 * y_ellipse_fine / span) ** 2, 0.0))
ax1.plot(y_ellipse_fine, ellipse_fine, "k--", linewidth=1.0, label="Elliptic (matched $c_\\mathrm{root}$)")
ax1.plot(y_full, chord_full, "o-", color="tab:blue", markersize=4, linewidth=1.5, label="Optimized")
ax1.set_xlabel("Spanwise station $y$ [m]")
ax1.set_ylabel("Chord [m]")
ax1.set_title(f"Chord Distribution\n{subtitle}")
ax1.legend(loc="lower center")
ax1.set_xlim(-span / 2 * 1.05, span / 2 * 1.05)
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
y_ell_fine = np.linspace(-span / 2, span / 2, 200)
c_ell_fine = c_root * np.sqrt(np.maximum(1.0 - (2.0 * y_ell_fine / span) ** 2, 0.0))
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

plt.show()
print("Plots saved to openmdao_vlm_optimization_out/")
