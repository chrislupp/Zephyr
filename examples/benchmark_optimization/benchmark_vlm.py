# VLM Philote chord optimization benchmark.
#
# Same problem as openmdao_vlm_optimization.py but with timing.
#
# Prerequisites:
#   1. Start the mesh server:    ./mesh_philote_server localhost:50052
#   2. Start the VLM server:     ./vlm_philote_server localhost:50051

import time
import numpy as np
import grpc
import openmdao.api as om
import philote_mdo.openmdao as pmom


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


span = 10.0
chord_init = 1.0
alpha_deg = 5.0
CL_target = 0.5
nx, ny = 6, 12
n_stations = ny + 1
n_cp = 5

yle = np.linspace(0.0, span / 2.0, n_stations)
ref_area = span * chord_init

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

print(f"\n{'='*50}")
print(f"VLM Philote Optimization Results")
print(f"{'='*50}")
print(f"  CL = {prob['CL'][0]:.6f}  (target: {CL_target})")
print(f"  CD = {prob['CD'][0]:.6f}")
print(f"  L/D = {prob['CL'][0] / prob['CD'][0]:.2f}")
print(f"  alpha = {np.rad2deg(prob['alpha'][0]):.3f} deg")
print(f"  planform area = {prob['reference_area'][0]:.4f} m^2")
print(f"\n  Iterations: {n_iter}")
print(f"  Total time: {t_total:.3f} s")
print(f"  Time/iteration: {t_total / n_iter * 1000:.1f} ms")
print(f"{'='*50}\n")
