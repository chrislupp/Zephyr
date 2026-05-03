# OpenAeroStruct chord optimization benchmark.
#
# Same problem as openmdao_vlm_optimization.py:
#   - Rectangular wing, span=10m, initial chord=1m
#   - nx=6 chordwise, ny=12 spanwise panels (half-span, symmetry)
#   - Minimize CD subject to CL=0.5, fixed planform area
#   - 5 chord control points + alpha as design variables
#   - SLSQP, tol=1e-8, maxiter=200

import time
import numpy as np
import openmdao.api as om
from openaerostruct.geometry.geometry_group import Geometry
from openaerostruct.aerodynamics.aero_groups import AeroPoint

span = 10.0
chord_init = 1.0
alpha_deg = 5.0
v = 25.0
rho = 1.225
CL_target = 0.5
ref_area = span * chord_init

nx = 6
ny = 12

# Generate half-span rectangular mesh (OAS convention: tip at y<0, root at y=0)
num_x = nx + 1
num_y = ny + 1
mesh = np.zeros((num_x, num_y, 3))
y_stations = np.linspace(-span / 2.0, 0.0, num_y)
for i in range(num_x):
    x_frac = i / nx
    mesh[i, :, 0] = x_frac * chord_init
    mesh[i, :, 1] = y_stations

surface = {
    "name": "wing",
    "symmetry": True,
    "S_ref_type": "projected",
    "mesh": mesh,
    "chord_cp": np.ones(5),
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
indep.add_output("v", val=v, units="m/s")
indep.add_output("alpha", val=alpha_deg, units="deg")
indep.add_output("Mach_number", val=0.0)
indep.add_output("re", val=1e6, units="1/m")
indep.add_output("rho", val=rho, units="kg/m**3")
indep.add_output("cg", val=np.zeros(3), units="m")

# Follow OAS walkthrough pattern exactly
geom_group = Geometry(surface=surface)
model.add_subsystem(surface["name"], geom_group)

aero_group = AeroPoint(surfaces=[surface])
point_name = "aero_point_0"
model.add_subsystem(point_name, aero_group,
                    promotes_inputs=["v", "alpha", "Mach_number", "re", "rho", "cg"])

name = surface["name"]
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

print(f"\n{'='*50}")
print(f"OpenAeroStruct Optimization Results")
print(f"{'='*50}")
print(f"  CL = {prob[point_name + '.CL'][0]:.6f}  (target: {CL_target})")
print(f"  CD = {prob[point_name + '.CD'][0]:.6f}")
print(f"  L/D = {prob[point_name + '.CL'][0] / prob[point_name + '.CD'][0]:.2f}")
print(f"  alpha = {prob['alpha'][0]:.3f} deg")
print(f"  S_ref = {prob[point_name + '.' + name + '.S_ref'][0]:.4f} m^2")
print(f"\n  Iterations: {n_iter}")
print(f"  Total time: {t_total:.3f} s")
print(f"  Time/iteration: {t_total / n_iter * 1000:.1f} ms")
print(f"{'='*50}\n")
