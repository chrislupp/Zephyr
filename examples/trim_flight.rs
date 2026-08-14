use pegasus_vlm::core::Vec3;
use pegasus_vlm::geometry::{Aircraft, Section, Surface};
use pegasus_vlm::solver::VortexLatticeSolver;
use pegasus_vlm::trim::{TrimPresets, TrimSolver};

fn main() -> pegasus_vlm::Result<()> {
    let mut wing = Surface::new("Wing");
    wing.set_chordwise_panel_count(5);
    wing.set_spanwise_panel_count(16);
    wing.add_section(Section::new(Vec3::new(0.0, -5.0, 0.0), 1.5, 0.0));
    wing.add_section(Section::new(Vec3::new(0.0, 5.0, 0.0), 1.5, 0.0));

    let mut tail = Surface::new("Horizontal tail");
    tail.set_chordwise_panel_count(3);
    tail.set_spanwise_panel_count(8);
    tail.add_section(Section::new(Vec3::new(4.0, -1.5, 0.0), 0.6, 0.0));
    tail.add_section(Section::new(Vec3::new(4.0, 1.5, 0.0), 0.4, 0.0));

    let mut aircraft = Aircraft::new("Wing-tail trim example");
    aircraft.set_reference_area(15.0);
    aircraft.set_reference_chord(1.5);
    aircraft.set_reference_span(10.0);
    aircraft.set_moment_reference_point(Vec3::new(0.375, 0.0, 0.0));
    aircraft.add_surface(wing);
    aircraft.add_surface(tail);

    let mut solver = VortexLatticeSolver::new();
    solver.set_aircraft(&aircraft)?;
    let result = TrimSolver::new().solve(
        &aircraft,
        &mut solver,
        &TrimPresets::level_flight(0.5),
        None,
    );
    println!("status: {:?}", result.status);
    println!("alpha: {:.4} deg", result.state.alpha.to_degrees());
    println!("CL: {:.6}", result.aerodynamic_results.coefficients.lift);
    println!("residual: {:.3e}", result.residual_norm);
    Ok(())
}
