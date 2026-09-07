//! Ported from the frozen v0.2.0 `mesh_node_gradient_test.cpp`.
//!
//! Adjoint gradients of the aerodynamic coefficients with respect to every
//! node of a directly supplied mesh, validated against re-solves.
//!
//! The frozen cases allowed a pass *rate* — 80% of components within tolerance
//! for the full-span mesh, 95% for the duplicated one — rather than requiring
//! every component to agree. The Rust adjoint does not need that concession:
//! every component agrees, with a worst relative error near 6e-7, which is the
//! truncation floor of the difference itself. These assert full agreement and
//! bound the worst error, and record the oracle's weaker bars so the
//! strengthening is visible rather than silent.

use pegasus_vlm::aero::ForceCalculator;
use pegasus_vlm::core::{DEGREES_TO_RADIANS, Vec3};
use pegasus_vlm::geometry::{MeshInput, SurfaceMesh};
use pegasus_vlm::gradients::GradientCalculator;
use pegasus_vlm::solver::{FlightCondition, VortexLatticeSolver};

/// A flat rectangular mesh, optionally mirrored about `y = 0`.
fn mesh_input(nx: usize, ny: usize, chord: f64, span: f64, duplicated: bool) -> MeshInput {
    let mut surface = SurfaceMesh::new("Wing", nx, ny);
    for i in 0..=nx {
        for j in 0..=ny {
            let x = chord * i as f64 / nx as f64;
            let y = span * j as f64 / ny as f64;
            surface.set_node(i, j, Vec3::new(x, y, 0.0)).unwrap();
        }
    }
    if duplicated {
        surface.y_duplicated = true;
        surface.y_duplicate_plane = 0.0;
    }
    let mut input = MeshInput::default();
    input.surfaces.push(surface);
    input.reference.area = chord * span;
    input.reference.chord = chord;
    input.reference.span = span;
    input.reference.moment_reference_point = Vec3::new(0.25 * chord, 0.5 * span, 0.0);
    input
}

fn condition() -> FlightCondition {
    FlightCondition {
        alpha: 5.0 * DEGREES_TO_RADIANS,
        ..FlightCondition::default()
    }
}

/// Lift and drag from a full re-solve of the given mesh.
fn coefficients(input: &MeshInput) -> (f64, f64) {
    let mut solver = VortexLatticeSolver::new();
    solver.set_mesh(input).unwrap();
    let solution = solver.solve(&condition()).unwrap();
    let results = ForceCalculator::new().compute_with_references(
        &input.reference,
        solver.mesh(),
        &solution,
        &condition(),
    );
    (results.coefficients.lift, results.coefficients.drag)
}

/// Differences every node coordinate and reports the fraction of components
/// that agree with the adjoint, plus the worst relative error seen.
fn validate(input: &MeshInput, nx: usize, ny: usize) -> (f64, f64) {
    let mut calculator = GradientCalculator::new();
    calculator.set_mesh(input).unwrap();
    let gradients = calculator
        .compute_mesh_node_gradients(&condition())
        .unwrap();
    assert_eq!(gradients.len(), 1);
    let lift = &gradients[0].lift;
    let drag = &gradients[0].drag;

    let step = 1.0e-6;
    let node_count = (nx + 1) * (ny + 1);
    let mut checked = 0usize;
    let mut agreed = 0usize;
    let mut worst = 0.0f64;

    for node in 0..node_count {
        for direction in 0..3 {
            let mut plus = input.clone();
            let mut minus = input.clone();
            plus.surfaces[0].nodes[node][direction] += step;
            minus.surfaces[0].nodes[node][direction] -= step;
            let (lift_plus, drag_plus) = coefficients(&plus);
            let (lift_minus, drag_minus) = coefficients(&minus);

            for (exact, forward, backward) in [
                (lift.nodes[node][direction], lift_plus, lift_minus),
                (drag.nodes[node][direction], drag_plus, drag_minus),
            ] {
                let difference = (forward - backward) / (2.0 * step);
                let tolerance = 1.0e-4 * (difference.abs() + 1.0e-6);
                let error = (exact - difference).abs();
                checked += 1;
                if error <= tolerance {
                    agreed += 1;
                }
                let scale = difference.abs().max(1.0e-6);
                worst = worst.max(error / scale);
            }
        }
    }
    (agreed as f64 / checked as f64, worst)
}

/// Ports `OutputStructure`.
#[test]
fn node_gradients_have_one_entry_per_node_and_surface() {
    let (nx, ny) = (2, 4);
    let input = mesh_input(nx, ny, 1.0, 4.0, false);
    let mut calculator = GradientCalculator::new();
    calculator.set_mesh(&input).unwrap();
    let gradients = calculator
        .compute_mesh_node_gradients(&condition())
        .unwrap();

    assert_eq!(gradients.len(), 1);
    assert_eq!(gradients[0].lift.nx, nx);
    assert_eq!(gradients[0].lift.ny, ny);
    assert_eq!(gradients[0].lift.nodes.len(), (nx + 1) * (ny + 1));
    assert_eq!(gradients[0].drag.nodes.len(), (nx + 1) * (ny + 1));
}

/// Ports `CLGradientFiniteDifference`.
///
/// The frozen bar was an 80% pass rate over all node components for both lift
/// and drag; every component agrees here.
#[test]
fn node_gradients_agree_with_resolves_on_a_full_span_mesh() {
    let (nx, ny) = (2, 3);
    let input = mesh_input(nx, ny, 1.0, 3.0, false);
    let (rate, worst) = validate(&input, nx, ny);
    // The oracle asked only for 0.8 here.
    assert_eq!(rate, 1.0, "worst relative error {worst:e}");
    assert!(worst < 1.0e-5, "worst relative error {worst:e}");
}

/// Ports `CLGradientYDuplicated`.
///
/// A mirrored mesh must do better, not worse: the frozen bar was 95%. Every node
/// contributes through both its own panel and its reflection, so a mishandled
/// mirror shows up as a large drop here while leaving the full-span case alone.
#[test]
fn node_gradients_agree_with_resolves_on_a_duplicated_mesh() {
    let (nx, ny) = (2, 3);
    let input = mesh_input(nx, ny, 1.0, 3.0, true);
    let (rate, worst) = validate(&input, nx, ny);
    // The oracle asked only for 0.95 here.
    assert_eq!(rate, 1.0, "worst relative error {worst:e}");
    assert!(worst < 1.0e-5, "worst relative error {worst:e}");
}

/// Ports `FullSpanVsYDuplicated`.
///
/// A mesh spanning `[-b/2, b/2]` outright and a half mesh mirrored about
/// `y = 0` describe the same wing, so their coefficients must agree. This is
/// what catches a duplication that builds the wrong geometry rather than the
/// wrong gradient.
#[test]
fn a_duplicated_half_mesh_matches_the_equivalent_full_span_mesh() {
    let (nx, ny) = (2, 3);
    let half_span = 1.5;

    let duplicated = mesh_input(nx, ny, 1.0, half_span, true);

    // The same wing meshed outright, spanning both halves.
    let mut surface = SurfaceMesh::new("Wing", nx, 2 * ny);
    for i in 0..=nx {
        for j in 0..=(2 * ny) {
            let x = i as f64 / nx as f64;
            let y = -half_span + 2.0 * half_span * j as f64 / (2 * ny) as f64;
            surface.set_node(i, j, Vec3::new(x, y, 0.0)).unwrap();
        }
    }
    let mut full = MeshInput::default();
    full.surfaces.push(surface);
    full.reference = duplicated.reference.clone();

    let (duplicated_lift, _) = coefficients(&duplicated);
    let (full_lift, _) = coefficients(&full);
    // Panel distributions differ slightly at the centreline, so this is a
    // physical agreement rather than a bit-for-bit one.
    assert!(
        (duplicated_lift - full_lift).abs() < 0.05 * full_lift.abs().max(1.0e-6),
        "duplicated {duplicated_lift} vs full span {full_lift}"
    );
}

/// Ports `NormalGradientKernel`.
///
/// The panel-normal JVP in isolation: perturbing a bound endpoint rotates the
/// panel normal, and the closed form must match a difference of the normal
/// itself.
#[test]
fn the_panel_normal_jvp_matches_a_centered_difference() {
    use pegasus_vlm::gradients::panel_normal_jvp;

    let left = Vec3::new(0.25, 0.0, 0.0);
    let right = Vec3::new(0.25, 1.0, 0.0);
    let step = 1.0e-7;

    // Normal of a flat panel with no twist and no camber, built the same way
    // the JVP differentiates it.
    let normal_of = |left: Vec3, right: Vec3| -> Vec3 {
        let bound = right - left;
        let spanwise = Vec3::new(0.0, bound.y, bound.z).normalize();
        Vec3::new(1.0, 0.0, 0.0).cross(&spanwise).normalize()
    };

    for direction_index in 0..3 {
        let mut direction = Vec3::zeros();
        direction[direction_index] = 1.0;
        let exact = panel_normal_jvp(
            &left,
            &right,
            &Vec3::zeros(),
            &direction,
            0.0,
            0.0,
            0.0,
            (0.0, 1.0),
            (0.0, 0.0),
        );
        let mut plus = right;
        let mut minus = right;
        plus[direction_index] += step;
        minus[direction_index] -= step;
        let difference = (normal_of(left, plus) - normal_of(left, minus)) / (2.0 * step);
        // The frozen tolerance.
        assert!(
            (exact - difference).norm() < 1.0e-5,
            "component {direction_index}: {exact:?} vs {difference:?}"
        );
    }
}

/// Far-field induced drag from a full re-solve, for the Trefftz-plane path.
fn induced_drag(input: &MeshInput) -> f64 {
    let mut solver = VortexLatticeSolver::new();
    solver.set_mesh(input).unwrap();
    let solution = solver.solve(&condition()).unwrap();
    let results = ForceCalculator::new().compute_with_references(
        &input.reference,
        solver.mesh(),
        &solution,
        &condition(),
    );
    results.far_field_induced_drag
}

/// The induced-drag node gradient comes from the Trefftz sweep rather than the
/// near-field force sweep, so it needs its own check against re-solves.
#[test]
fn induced_drag_node_gradients_agree_with_resolves() {
    let (nx, ny) = (2, 3);
    let input = mesh_input(nx, ny, 1.0, 3.0, false);
    let mut calculator = GradientCalculator::new();
    calculator.set_mesh(&input).unwrap();
    let gradients = calculator
        .compute_mesh_node_gradients(&condition())
        .unwrap();
    let exact = &gradients[0].induced_drag;

    let step = 1.0e-6;
    let mut worst: f64 = 0.0;
    for node in 0..(nx + 1) * (ny + 1) {
        for direction in 0..3 {
            let mut plus = input.clone();
            let mut minus = input.clone();
            plus.surfaces[0].nodes[node][direction] += step;
            minus.surfaces[0].nodes[node][direction] -= step;
            let difference = (induced_drag(&plus) - induced_drag(&minus)) / (2.0 * step);
            let error = (exact.nodes[node][direction] - difference).abs();
            worst = worst.max(error / difference.abs().max(1.0e-6));
        }
    }
    assert!(worst < 1.0e-5, "worst relative error {worst}");
}

/// Same, with the mirrored half folded back onto the original surface.
#[test]
fn induced_drag_node_gradients_agree_on_a_duplicated_mesh() {
    let (nx, ny) = (2, 3);
    let input = mesh_input(nx, ny, 1.0, 3.0, true);
    let mut calculator = GradientCalculator::new();
    calculator.set_mesh(&input).unwrap();
    let gradients = calculator
        .compute_mesh_node_gradients(&condition())
        .unwrap();
    let exact = &gradients[0].induced_drag;

    let step = 1.0e-6;
    let mut worst: f64 = 0.0;
    for node in 0..(nx + 1) * (ny + 1) {
        for direction in 0..3 {
            let mut plus = input.clone();
            let mut minus = input.clone();
            plus.surfaces[0].nodes[node][direction] += step;
            minus.surfaces[0].nodes[node][direction] -= step;
            let difference = (induced_drag(&plus) - induced_drag(&minus)) / (2.0 * step);
            let error = (exact.nodes[node][direction] - difference).abs();
            worst = worst.max(error / difference.abs().max(1.0e-6));
        }
    }
    assert!(worst < 1.0e-5, "worst relative error {worst}");
}
