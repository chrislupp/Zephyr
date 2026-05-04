# Mesh scaling benchmark: VLM Philote vs OpenAeroStruct.
#
# Sweeps spanwise panel count (ny) with fixed nx=6 and plots:
#   1. Wall time per iteration vs total panels
#   2. Total optimization wall time vs total panels
#
# Prerequisites:
#   1. Start the mesh server:    ./mesh_philote_server localhost:50052
#   2. Start the VLM server:     ./vlm_philote_server localhost:50051

import time
import numpy as np
import grpc
import openmdao.api as om
import philote_mdo.openmdao as pmom
from openaerostruct.geometry.geometry_group import Geometry
from openaerostruct.aerodynamics.aero_groups import AeroPoint

# ── Shared parameters ────────────────────────────────────────────────────

span = 10.0
chord_init = 1.0
alpha_deg = 5.0
CL_target = 0.5
ref_area = span * chord_init
nx = 6
n_cp = 5
ny_values = [8, 12, 16, 24, 32, 48]


# ── Helper components for VLM Philote ────────────────────────────────────

class AkimaChord(om.ExplicitComponent):
    def initialize(self):
        self.options.declare("n_cp", types=int)
        self.options.declare("n_stations", types=int)
        self.options.declare("y_stations", types=np.ndarray)

    def setup(self):
        self.add_input("chord_cp", shape=(self.options["n_cp"],), units="m")
        self.add_output("chord", shape=(self.options["n_stations"],), units="m")

    def setup_partials(self):
        self.declare_partials("chord", "chord_cp")

    def compute(self, inputs, outputs):
        from scipy.interpolate import Akima1DInterpolator
        n_cp = self.options["n_cp"]
        y_st = self.options["y_stations"]
        y_cp = np.linspace(0.0, y_st[-1], n_cp)
        outputs["chord"] = Akima1DInterpolator(y_cp, inputs["chord_cp"])(y_st)

    def compute_partials(self, inputs, partials):
        from scipy.interpolate import Akima1DInterpolator
        n_cp, n_st = self.options["n_cp"], self.options["n_stations"]
        y_st = self.options["y_stations"]
        y_cp = np.linspace(0.0, y_st[-1], n_cp)
        cp = inputs["chord_cp"].ravel()
        jac = np.zeros((n_st, n_cp))
        h = 1e-8
        for k in range(n_cp):
            cp_p, cp_m = cp.copy(), cp.copy()
            cp_p[k] += h; cp_m[k] -= h
            jac[:, k] = (Akima1DInterpolator(y_cp, cp_p)(y_st)
                         - Akima1DInterpolator(y_cp, cp_m)(y_st)) / (2*h)
        partials["chord", "chord_cp"] = jac


# ── VLM Philote benchmark ───────────────────────────────────────────────

def run_vlm(ny):
    n_stations = ny + 1
    yle = np.linspace(0.0, span / 2.0, n_stations)

    prob = om.Problem()
    model = prob.model

    model.add_subsystem("chord_interp",
        AkimaChord(n_cp=n_cp, n_stations=n_stations, y_stations=yle),
        promotes_inputs=[("chord_cp", "chord_cp")], promotes_outputs=["chord"])

    mesh_channel = grpc.insecure_channel("localhost:50052")
    model.add_subsystem("mesh_comp",
        pmom.RemoteExplicitComponent(channel=mesh_channel, nx=nx, ny=ny,
            cosine_chordwise_spacing="true", cosine_spanwise_spacing="false",
            y_duplicated="true"),
        promotes_inputs=["chord", "twist", "xle", "yle", "zle"],
        promotes_outputs=["reference_area", "reference_span", "reference_chord"])

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

    prob.driver = om.ScipyOptimizeDriver()
    prob.driver.options["optimizer"] = "SLSQP"
    prob.driver.options["tol"] = 1e-8
    prob.driver.options["maxiter"] = 200

    model.add_design_var("chord_cp", lower=0.1, upper=3.0)
    model.add_design_var("alpha", lower=np.deg2rad(-10.0), upper=np.deg2rad(15.0))
    model.add_objective("CD")
    model.add_constraint("CL", equals=CL_target)
    model.add_constraint("reference_area", equals=ref_area)

    prob.setup()
    prob["chord_cp"] = np.full(n_cp, chord_init)
    prob["twist"] = np.zeros(n_stations)
    prob["xle"] = np.zeros(n_stations)
    prob["yle"] = yle
    prob["zle"] = np.zeros(n_stations)
    # reference_area, reference_chord, reference_span computed by mesh_comp
    prob["moment_ref_x"] = 0.25 * chord_init
    prob["moment_ref_y"] = 0.0
    prob["moment_ref_z"] = 0.0
    prob["alpha"] = np.deg2rad(alpha_deg)
    prob["beta"] = 0.0
    prob["mach"] = 0.0
    prob["dynamic_pressure"] = 0.5 * 1.225 * 25.0**2
    prob["roll_rate"] = 0.0
    prob["pitch_rate"] = 0.0
    prob["yaw_rate"] = 0.0

    prob.run_model()
    t0 = time.perf_counter()
    prob.run_driver()
    t_total = time.perf_counter() - t0
    n_iter = prob.driver.iter_count
    return t_total, n_iter


# ── OAS benchmark ────────────────────────────────────────────────────────

def run_oas(ny):
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
        "chord_cp": np.ones(n_cp),
        "CL0": 0.0,
        "CD0": 0.0,
        "with_viscous": False,
        "with_wave": False,
        "k_lam": 0.05,
        "t_over_c_cp": np.array([0.12]),
        "c_max_t": 0.303,
    }

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

    prob.driver = om.ScipyOptimizeDriver()
    prob.driver.options["optimizer"] = "SLSQP"
    prob.driver.options["tol"] = 1e-8
    prob.driver.options["maxiter"] = 200

    model.add_design_var("wing.chord_cp", lower=0.1, upper=3.0)
    model.add_design_var("alpha", lower=-10.0, upper=15.0, units="deg")
    model.add_objective(point_name + ".CD")
    model.add_constraint(point_name + ".CL", equals=CL_target)
    model.add_constraint(point_name + "." + name + ".S_ref", equals=ref_area)

    prob.setup()
    prob.run_model()
    t0 = time.perf_counter()
    prob.run_driver()
    t_total = time.perf_counter() - t0
    n_iter = prob.driver.iter_count
    return t_total, n_iter


# ── Run sweep ────────────────────────────────────────────────────────────

vlm_results = []
oas_results = []

for ny in ny_values:
    n_panels = nx * 2 * ny  # full-span panel count
    print(f"\nny={ny}  ({n_panels} panels)")

    print("  VLM ...", end=" ", flush=True)
    t, n = run_vlm(ny)
    vlm_results.append((n_panels, t, n))
    print(f"{t:.2f}s, {n} iters, {t/n*1000:.1f} ms/iter")

    print("  OAS ...", end=" ", flush=True)
    t, n = run_oas(ny)
    oas_results.append((n_panels, t, n))
    print(f"{t:.2f}s, {n} iters, {t/n*1000:.1f} ms/iter")

vlm_results = np.array(vlm_results)
oas_results = np.array(oas_results)

# ── Plotting ─────────────────────────────────────────────────────────────

import matplotlib.pyplot as plt

panels_vlm = vlm_results[:, 0]
panels_oas = oas_results[:, 0]

# Figure 1: Wall time per iteration
fig1, ax1 = plt.subplots(figsize=(7, 5))
ax1.plot(panels_vlm, vlm_results[:, 1] / vlm_results[:, 2] * 1000,
         "o-", color="tab:blue", label="VLM Philote")
ax1.plot(panels_oas, oas_results[:, 1] / oas_results[:, 2] * 1000,
         "s-", color="tab:orange", label="OpenAeroStruct")
ax1.set_xlabel("Total panels (full-span)")
ax1.set_ylabel("Wall time per iteration [ms]")
ax1.set_title("Time per Iteration vs Mesh Size")
ax1.legend()
ax1.grid(True, alpha=0.3)
fig1.tight_layout()
fig1.savefig("benchmark_scaling_out/time_per_iteration.png", dpi=150)

# Figure 2: Total optimization time
fig2, ax2 = plt.subplots(figsize=(7, 5))
ax2.plot(panels_vlm, vlm_results[:, 1],
         "o-", color="tab:blue", label="VLM Philote")
ax2.plot(panels_oas, oas_results[:, 1],
         "s-", color="tab:orange", label="OpenAeroStruct")
ax2.set_xlabel("Total panels (full-span)")
ax2.set_ylabel("Total wall time [s]")
ax2.set_title("Optimization Wall Time vs Mesh Size")
ax2.legend()
ax2.grid(True, alpha=0.3)
ax2.set_xscale("log")
ax2.set_yscale("log")
fig2.tight_layout()
fig2.savefig("benchmark_scaling_out/total_time.png", dpi=150)

plt.show()
print("Plots saved to benchmark_scaling_out/")
