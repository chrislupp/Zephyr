#!/usr/bin/env python3
"""Gradient-based sizing with internal trim through explicit Philote RPCs.

Start a fresh philote_vlm_server and build philote_trim_optimization, then run
this script with NumPy and SciPy installed. The objective gradient includes the
converged trim response. Finite differences are used only for verification.
"""
import argparse
import functools
import json
from pathlib import Path
import subprocess

import numpy as np
import scipy
from scipy.optimize import minimize


class TrimmedModel:
    def __init__(self, binary, endpoint):
        self.process = subprocess.Popen(
            [str(binary), endpoint], stdin=subprocess.PIPE, stdout=subprocess.PIPE,
            text=True, bufsize=1,
        )
        self.evaluations = 0
        self.max_force_residual = 0.0
        self.max_moment_residual = 0.0

    @functools.lru_cache(maxsize=256)
    def _evaluate(self, design):
        self.process.stdin.write(json.dumps(design) + "\n")
        self.process.stdin.flush()
        line = self.process.stdout.readline()
        if not line:
            raise RuntimeError("Philote optimization bridge exited")
        result = json.loads(line)
        if "error" in result:
            raise RuntimeError(result["error"])
        self.evaluations += 1
        self.max_force_residual = max(self.max_force_residual, abs(result["force_residual"]))
        self.max_moment_residual = max(self.max_moment_residual, abs(result["moment_residual"]))
        return result

    def evaluate(self, design):
        return self._evaluate(tuple(float(value) for value in design))

    def close(self):
        self.process.stdin.close()
        code = self.process.wait(timeout=30)
        if code:
            raise RuntimeError(f"Philote optimization bridge exited with {code}")


def check_gradients(model, design):
    value = model.evaluate(design)
    analytical = np.vstack([value["objective_gradient"], value["alpha_gradient"], value["elevator_gradient"]])
    numerical = np.empty_like(analytical)
    for column in range(2):
        h = 1e-5 * max(1.0, abs(design[column]))
        plus = np.array(design, dtype=float)
        minus = plus.copy()
        plus[column] += h
        minus[column] -= h
        p, m = model.evaluate(plus), model.evaluate(minus)
        numerical[:, column] = (
            np.array([p["objective"], p["alpha"], p["elevator"]])
            - np.array([m["objective"], m["alpha"], m["elevator"]])
        ) / (2 * h)
    error = float(np.max(np.abs(analytical - numerical) / (1 + np.abs(numerical))))
    if error > 1e-5:
        raise RuntimeError(f"analytical gradient check failed: {error:.3e}")
    return {"normalized_max_error": error, "analytical": analytical.tolist(), "finite_difference": numerical.tolist()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    root = Path(__file__).resolve().parents[1]
    parser.add_argument("--binary", type=Path, default=root / "target/debug/examples/philote_trim_optimization")
    parser.add_argument("--endpoint", default="http://127.0.0.1:50051")
    parser.add_argument("--output", type=Path, default=Path("trim_optimization_results.json"))
    parser.add_argument("--initial", type=float, nargs=2, default=[4.0, 0.65], metavar=("HALFSPAN", "TAIL_CHORD"))
    args = parser.parse_args()
    model = TrimmedModel(args.binary, args.endpoint)
    bounds = [(2.2, 5.0), (0.3, 1.2)]
    initial = np.array(args.initial)
    objective_scale = 30.0
    history = []
    try:
        before = model.evaluate(initial)
        initial_check = check_gradients(model, initial)
        result = minimize(
            lambda x: model.evaluate(x)["objective"] / objective_scale, initial,
            jac=lambda x: np.array(model.evaluate(x)["objective_gradient"]) / objective_scale,
            method="SLSQP", bounds=bounds,
            callback=lambda x: history.append(model.evaluate(x)),
            options={"ftol": 1e-8, "maxiter": 100, "disp": True},
        )
        after = model.evaluate(result.x)
        final_check = check_gradients(model, result.x)
        if not result.success:
            raise RuntimeError(f"optimization did not converge: {result.message}")
        if after["objective"] >= before["objective"]:
            raise RuntimeError("optimization did not improve the objective")
        if model.max_force_residual > 1e-6 or model.max_moment_residual > 1e-6:
            raise RuntimeError("trim was not maintained across evaluations")
        report = {"scipy_version": scipy.__version__, "optimizer": "SLSQP", "success": bool(result.success),
                  "message": str(result.message), "iterations": int(result.nit), "evaluations": model.evaluations,
                  "initial": before, "final": after, "bounds": bounds,
                  "objective_improvement_percent": 100 * (before["objective"] - after["objective"]) / before["objective"],
                  "initial_gradient_check": initial_check, "final_gradient_check": final_check,
                  "maximum_force_residual": model.max_force_residual,
                  "maximum_moment_residual": model.max_moment_residual, "history": history}
        args.output.write_text(json.dumps(report, indent=2) + "\n")
        print(f"Objective: {before['objective']:.6f} -> {after['objective']:.6f} N")
        print(f"Design [halfspan, tail chord]: {initial} -> {result.x} m")
        print(f"Trim [alpha, elevator]: [{after['alpha']:.6f}, {after['elevator']:.6f}] rad")
        print(f"Maximum balance residuals: {model.max_force_residual:.3e} N, {model.max_moment_residual:.3e} N m")
        print(f"Gradient errors: {initial_check['normalized_max_error']:.3e}, {final_check['normalized_max_error']:.3e}")
        print(f"Report: {args.output}")
    finally:
        model.close()


if __name__ == "__main__":
    main()
