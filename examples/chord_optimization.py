#!/usr/bin/env python3
#
#     Project Zephyr
#
#     Copyright 2024-2026, Christopher A. Lupp
#
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     This work has been cleared for public release, distribution unlimited,
#     case number: AFRL-2026-2198.
#
#     The views expressed are those of the authors and do not reflect the
#     official guidance or position of the United States Government, the
#     Department of Defense or of the United States Air Force.
#
#     Statement from DoD: The Appearance of external hyperlinks does not
#     constitute endorsement by the United States Department of Defense (DoD)
#     of the linked websites, of the information, products, or services
#     contained therein. The DoD does not exercise any editorial, security, or
#     other control over the information you may find at these locations.
"""Minimum-induced-drag chord distribution, in OpenMDAO through Philote.

A rectangular wing is handed to IPOPT with its chord distribution as the design
variable and its far-field induced drag as the objective. It comes back
elliptical, which is the answer lifting-line theory gives for minimum induced
drag at fixed span and lift, so the run is a check on the analytical gradients
as much as it is an optimization.

    minimize    CDi
    with
    respect to  chord_cp[6]         B-spline control points of the chord
                alpha               angle of attack
    subject to  CL = 0.5
                reference_area = 10
                0.05 <= chord_cp <= 2.5

The chords are driven by a B-spline rather than being design variables one per
station. A chord per station leaves the optimizer free directions it should not
have: two neighbouring stations bound the same lattice strip, which sees only
their average, so a sawtooth in the station chords barely changes the loading
and the optimizer is content to leave one there. The result optimizes correctly
and looks terrible. Six control points span the same shapes without the
sawtooth, and are what a wing parameterization looks like in practice.

Two Zephyr disciplines are used, each behind its own gRPC server:

    ZephyrMesh  chord, twist, xle, zle, eta, halfspan -> mesh, reference_*
    ZephyrVlm   mesh, reference_*, flight condition   -> CL, CDi

Both report analytical partials, so OpenMDAO assembles the total derivatives
from the chain rule and nothing is finite-differenced. Start the servers first:

    cargo run --release --features philote --example philote_mesh_server
    cargo run --release --features philote --example philote_vlm_server

then run this script.
"""

import argparse
import json
import re
from pathlib import Path

import grpc
import numpy as np
import openmdao.api as om
import philote_mdo.openmdao as pmom

# Span and reference area are held fixed: induced drag at fixed lift is a
# question about span, and the area constraint is what keeps the problem a
# question about shape rather than about scale.
SPAN = 10.0
AREA = 10.0
ASPECT_RATIO = SPAN**2 / AREA
TARGET_LIFT = 0.5

# Chordwise panels, and spanwise panels per semispan. The mesh discipline takes
# one chord per spanwise node, so the lattice fixes the number of stations the
# spline is evaluated at.
NX = 4
NY = 20
STATIONS = NY + 1

# B-spline control points driving those stations. Enough to reach an ellipse,
# too few to hold a sawtooth: six is the fewest that gives a strictly monotone
# planform here, and adding more changes the drag in the eighth decimal.
CONTROL_POINTS = 6
SPLINE_ORDER = 4

# A B-spline curve stays inside the convex hull of its control points, so
# bounding the control points bounds the chord.
CHORD_BOUNDS = (0.05, 2.5)


def elliptical_chords(eta):
    """The elliptical planform of the same span and area."""
    root = 4.0 * AREA / (np.pi * SPAN)
    return root * np.sqrt(np.maximum(0.0, 1.0 - eta**2))


def build_problem(mesh_address, vlm_address):
    """The OpenMDAO problem: two remote disciplines and one driver."""
    model = om.Group()

    # Design variables and the planform parameters held fixed.
    independents = om.IndepVarComp()
    independents.add_output("chord_cp", np.full(CONTROL_POINTS, AREA / SPAN), units="m")
    independents.add_output("twist", np.zeros(STATIONS), units="rad")
    independents.add_output("zle", np.zeros(STATIONS), units="m")
    independents.add_output("eta", np.linspace(0.0, 1.0, STATIONS))
    independents.add_output("halfspan", SPAN / 2.0, units="m")
    independents.add_output("alpha", np.radians(5.0), units="rad")
    model.add_subsystem("design", independents, promotes=["*"])

    # The design variables are control points; the mesh discipline wants a chord
    # at every spanwise station. The spline is one more link in the chain
    # OpenMDAO differentiates.
    spline = om.SplineComp(
        method="bsplines",
        x_interp_val=np.linspace(0.0, 1.0, STATIONS),
        num_cp=CONTROL_POINTS,
        interp_options={"order": SPLINE_ORDER},
    )
    spline.add_spline(
        y_cp_name="chord_cp",
        y_interp_name="chord_interp",
        y_cp_val=np.full(CONTROL_POINTS, AREA / SPAN),
        y_units="m",
    )
    # SplineComp carries a leading "one spline" dimension, so its output is
    # (1, STATIONS) where everything downstream wants (STATIONS,). Say so,
    # rather than leaving OpenMDAO to match them by size.
    model.add_subsystem("chord_spline", spline, promotes_inputs=["chord_cp"])
    model.connect("chord_spline.chord_interp", "chord", src_indices=om.slicer[0, :])

    # Hold the quarter-chord line straight and unswept as the chords change.
    # This is the geometry the elliptical result is derived for.
    model.add_subsystem(
        "quarter_chord",
        om.ExecComp(
            "xle = -0.25 * chord",
            xle={"shape": STATIONS, "units": "m"},
            chord={"shape": STATIONS, "units": "m"},
        ),
        promotes=["*"],
    )

    # Zephyr's mesh discipline: section parameters to a structured mesh and the
    # reference quantities that normalize the coefficients.
    model.add_subsystem(
        "geometry",
        pmom.RemoteExplicitComponent(
            channel=grpc.insecure_channel(mesh_address),
            nx=NX,
            ny=NY,
            y_duplicated=True,
        ),
        promotes=["*"],
    )

    # Zephyr's VLM discipline. The selective-analysis options ask for the two
    # outputs this problem needs and for exactly the partials OpenMDAO will
    # use; everything else is never computed. "Wing" is the surface's node
    # input, so it is also the name of its mesh-node partials.
    model.add_subsystem(
        "aerodynamics",
        pmom.RemoteExplicitComponent(
            channel=grpc.insecure_channel(vlm_address),
            num_surfaces=1,
            surface_names="Wing",
            surface_nx=str(NX),
            surface_ny=str(NY),
            surface_y_duplicated="true",
            analysis_outputs="CL,CDi",
            analysis_gradients={
                "CL": ["alpha", "reference_area", "Wing"],
                "CDi": ["alpha", "reference_area", "Wing"],
            },
        ),
        promotes=["*"],
    )

    # reference_area, reference_chord and reference_span are promoted out of
    # the mesh discipline and into the VLM discipline, so they connect by name.
    # The mesh itself does not: name it explicitly.
    model.connect("mesh", "Wing")

    problem = om.Problem(model=model)
    problem.driver = om.pyOptSparseDriver(optimizer="IPOPT")
    problem.driver.options["print_results"] = False
    problem.driver.opt_settings.update(
        {
            "tol": 1e-8,
            "max_iter": 200,
            # Only first derivatives are available, so IPOPT builds curvature
            # itself rather than asking for a Hessian.
            "hessian_approximation": "limited-memory",
            "mu_strategy": "adaptive",
            "print_level": 0,
            "sb": "yes",
        }
    )

    problem.model.add_design_var(
        "chord_cp", lower=CHORD_BOUNDS[0], upper=CHORD_BOUNDS[1], units="m"
    )
    problem.model.add_design_var("alpha", lower=-0.3, upper=0.3, ref=0.1, units="rad")
    problem.model.add_objective("CDi", ref=0.01)
    problem.model.add_constraint("CL", equals=TARGET_LIFT, ref=TARGET_LIFT)
    problem.model.add_constraint("reference_area", equals=AREA, ref=AREA, units="m**2")
    return problem


def set_flight_condition(problem):
    """Every scalar the VLM discipline declares has to be given a value.

    OpenMDAO defaults an unconnected input to one, and a Mach number of one is
    rejected by the solver rather than silently clamped, so none of these can
    be left to the default.
    """
    problem.set_val("beta", 0.0, units="rad")
    problem.set_val("mach", 0.0)
    problem.set_val("dynamic_pressure", 1.0, units="Pa")
    problem.set_val("roll_rate", 0.0)
    problem.set_val("pitch_rate", 0.0)
    problem.set_val("yaw_rate", 0.0)
    problem.set_val("moment_ref_x", 0.0, units="m")
    problem.set_val("moment_ref_y", 0.0, units="m")
    problem.set_val("moment_ref_z", 0.0, units="m")


def trim_alpha(problem):
    """Put the wing on the lift constraint by solving for alpha.

    The baseline has to satisfy `CL = 0.5` before its drag can be compared with
    the optimum's, because induced drag goes as the square of lift. Newton on
    `dCL/dalpha`, which `compute_totals` gets from the same analytical chain the
    driver will use, converges in one step for a flat untwisted wing.
    """
    for _ in range(10):
        problem.run_model()
        residual = problem.get_val("CL")[0] - TARGET_LIFT
        if abs(residual) < 1e-12:
            break
        slope = problem.compute_totals(of=["CL"], wrt=["alpha"])["CL", "alpha"][0, 0]
        alpha = problem.get_val("alpha", units="rad")[0]
        problem.set_val("alpha", alpha - residual / slope, units="rad")


def snapshot(problem):
    """The quantities worth recording at a design point."""
    return {
        "induced_drag": float(problem.get_val("CDi")[0]),
        "lift": float(problem.get_val("CL")[0]),
        "alpha_degrees": float(np.degrees(problem.get_val("alpha", units="rad")[0])),
        "reference_area": float(problem.get_val("reference_area", units="m**2")[0]),
        "chords": problem.get_val("chord", units="m").tolist(),
        "control_points": problem.get_val("chord_cp", units="m").tolist(),
    }


# One row of IPOPT's iteration table: iteration number (with an "r" during
# feasibility restoration), objective, primal infeasibility, dual infeasibility.
IPOPT_ROW = re.compile(r"^\s*(\d+)r?\s+(\S+)\s+(\S+)\s+(\S+)\s")


def read_ipopt_iterations(path):
    """IPOPT's own convergence measures, read from its output file.

    ``inf_pr`` is the primal infeasibility, the largest constraint violation:
    feasibility. ``inf_du`` is the dual infeasibility, the largest component
    of the gradient of the Lagrangian: optimality. Both are in the scaled
    quantities IPOPT works with.
    """
    iterations = []
    for line in Path(path).read_text().splitlines():
        match = IPOPT_ROW.match(line)
        if not match:
            continue
        try:
            row = [float(value) for value in match.groups()]
        except ValueError:
            continue
        iterations.append(
            {
                "iteration": int(row[0]),
                "objective": row[1],
                "feasibility": row[2],
                "optimality": row[3],
            }
        )
    return iterations


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh-server", default="localhost:50052")
    parser.add_argument("--vlm-server", default="localhost:50051")
    parser.add_argument("--output", type=Path, default=Path("chord_optimization_results.json"))
    parser.add_argument(
        "--check-partials",
        action="store_true",
        help="compare the total derivatives against finite differences and stop",
    )
    args = parser.parse_args()

    problem = build_problem(args.mesh_server, args.vlm_server)
    problem.setup(mode="rev")
    set_flight_condition(problem)

    # Start from a rectangular wing that already carries the target lift, so the
    # drag reduction reported at the end is a like-for-like comparison.
    trim_alpha(problem)
    initial = snapshot(problem)
    eta = problem.get_val("eta")
    reference = elliptical_chords(eta)

    if args.check_partials:
        # The whole point of the exercise: the chain through two remote
        # disciplines, against finite differences across the same chain.
        problem.check_totals(
            of=["CDi", "CL", "reference_area"],
            wrt=["chord_cp", "alpha"],
            method="fd",
            step=1e-6,
            compact_print=True,
        )
        return

    failed = problem.run_driver()
    problem.run_model()
    final = snapshot(problem)
    # OpenMDAO directs IPOPT's iteration table to IPOPT.out in the problem's
    # outputs directory; the convergence history is read back from there.
    iterations = read_ipopt_iterations(Path(problem.get_outputs_dir()) / "IPOPT.out")

    chords = np.asarray(final["chords"])
    print()
    print(f"optimization {'FAILED' if failed else 'succeeded'} "
          f"in {iterations[-1]['iteration']} IPOPT iterations")
    print(f"final feasibility {iterations[-1]['feasibility']:.1e}, "
          f"optimality {iterations[-1]['optimality']:.1e}")
    print(f"CDi   {initial['induced_drag']:.8f} -> {final['induced_drag']:.8f}"
          f"   ({100.0 * (final['induced_drag'] / initial['induced_drag'] - 1.0):+.2f} %)")
    print(f"alpha {initial['alpha_degrees']:.3f} deg -> {final['alpha_degrees']:.3f} deg")
    print(f"CL    {final['lift']:.10f} (target {TARGET_LIFT})")
    print(f"S     {final['reference_area']:.10f} (target {AREA})")
    print(f"ideal CDi = CL^2 / (pi AR) = {TARGET_LIFT**2 / (np.pi * ASPECT_RATIO):.8f}")
    print()
    print("   y   initial  optimized  elliptical")
    for station in range(STATIONS):
        print(f"{eta[station] * SPAN / 2:5.2f}  {initial['chords'][station]:8.4f}"
              f"  {chords[station]:9.4f}  {reference[station]:10.4f}")

    report = {
        "problem": {
            "span": SPAN,
            "area": AREA,
            "aspect_ratio": ASPECT_RATIO,
            "target_lift": TARGET_LIFT,
            "stations": STATIONS,
            "control_points": CONTROL_POINTS,
            "spline_order": SPLINE_ORDER,
            "chord_bounds": list(CHORD_BOUNDS),
            "spanwise_panels": NY,
            "chordwise_panels": NX,
        },
        "succeeded": not failed,
        "stations": (eta * SPAN / 2).tolist(),
        "initial_chords": initial["chords"],
        "optimized_chords": chords.tolist(),
        "elliptical_chords": reference.tolist(),
        "initial": initial,
        "optimized": final,
        "ideal_induced_drag": TARGET_LIFT**2 / (np.pi * ASPECT_RATIO),
        "ipopt_iterations": iterations,
    }
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(f"\nwrote {args.output}")


if __name__ == "__main__":
    main()
