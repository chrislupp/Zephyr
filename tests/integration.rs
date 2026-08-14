use approx::assert_abs_diff_eq;
use pegasus_vlm::aero::ForceCalculator;
use pegasus_vlm::core::{DEGREES_TO_RADIANS, Vec3};
use pegasus_vlm::geometry::{Aircraft, Section, Surface};
use pegasus_vlm::gradients::{AeroOutput, TotalDerivative};
use pegasus_vlm::solver::{FlightCondition, VortexLatticeSolver};

fn rectangular_wing() -> Aircraft {
    let mut wing = Surface::new("Wing");
    wing.set_chordwise_panel_count(3);
    wing.set_spanwise_panel_count(8);
    wing.add_section(Section::new(Vec3::new(0.0, -2.0, 0.0), 1.0, 0.0));
    wing.add_section(Section::new(Vec3::new(0.0, 2.0, 0.0), 1.0, 0.0));
    let mut aircraft = Aircraft::new("Rectangular wing");
    aircraft.set_reference_area(4.0);
    aircraft.set_reference_chord(1.0);
    aircraft.set_reference_span(4.0);
    aircraft.set_moment_reference_point(Vec3::new(0.25, 0.0, 0.0));
    aircraft.add_surface(wing);
    aircraft
}

fn coefficient(aircraft: &Aircraft, condition: &FlightCondition, output: AeroOutput) -> f64 {
    let mut solver = VortexLatticeSolver::new();
    solver.set_aircraft(aircraft).unwrap();
    let solution = solver.solve(condition).unwrap();
    let result =
        ForceCalculator::new().compute_for_aircraft(aircraft, solver.mesh(), &solution, condition);
    match output {
        AeroOutput::Lift => result.coefficients.lift,
        AeroOutput::Drag => result.coefficients.drag,
        AeroOutput::Side => result.coefficients.side,
        AeroOutput::Roll => result.coefficients.roll,
        AeroOutput::Pitch => result.coefficients.pitch,
        AeroOutput::Yaw => result.coefficients.yaw,
        AeroOutput::InducedDrag => result.far_field_induced_drag,
    }
}

#[test]
fn rectangular_wing_solves_and_generates_positive_lift() {
    let aircraft = rectangular_wing();
    let condition = FlightCondition {
        alpha: 5.0 * DEGREES_TO_RADIANS,
        ..FlightCondition::default()
    };
    let mut solver = VortexLatticeSolver::new();
    solver.set_aircraft(&aircraft).unwrap();
    let solution = solver.solve(&condition).unwrap();
    assert!(solution.is_valid);
    assert_eq!(solution.circulation.len(), 24);
    let result = ForceCalculator::new().compute_for_aircraft(
        &aircraft,
        solver.mesh(),
        &solution,
        &condition,
    );
    assert!(result.is_valid);
    assert!(result.coefficients.lift > 0.0);
}

#[test]
fn adjoint_alpha_matches_full_resolves_for_all_primary_coefficients() {
    let aircraft = rectangular_wing();
    let condition = FlightCondition {
        alpha: 4.0 * DEGREES_TO_RADIANS,
        beta: 2.0 * DEGREES_TO_RADIANS,
        ..FlightCondition::default()
    };
    let mut solver = VortexLatticeSolver::new();
    solver.set_aircraft(&aircraft).unwrap();
    solver.solve(&condition).unwrap();
    let step = 1.0e-6;
    for output in [
        AeroOutput::Lift,
        AeroOutput::Drag,
        AeroOutput::Side,
        AeroOutput::Roll,
        AeroOutput::Pitch,
        AeroOutput::Yaw,
    ] {
        let analytical = TotalDerivative::compute(&solver, &condition, output)
            .unwrap()
            .by_alpha;
        let mut plus = condition.clone();
        plus.alpha += step;
        let mut minus = condition.clone();
        minus.alpha -= step;
        let finite_difference = (coefficient(&aircraft, &plus, output)
            - coefficient(&aircraft, &minus, output))
            / (2.0 * step);
        assert_abs_diff_eq!(analytical, finite_difference, epsilon = 2.0e-5);
    }
}
