//! Reverse-mode strip-load sensitivities to the mesh nodes.
//!
//! `StripLoadAdjoint` exists for aeroelastic coupling, where the structure needs
//! the derivative of each strip's resultant force and moment — not of an
//! aircraft-level coefficient — with respect to a mesh that moves with the
//! deformation. These difference the same functional the adjoint contracts,
//! re-solving the lattice at each perturbed mesh, so the check covers the
//! implicit circulation path as well as the explicit Kutta–Joukowski one.

use pegasus_vlm::aero::ForceCalculator;
use pegasus_vlm::core::{DEGREES_TO_RADIANS, Vec3};
use pegasus_vlm::geometry::{MeshInput, SurfaceMesh};
use pegasus_vlm::gradients::{StripLoadAdjoint, StripLoadCotangent};
use pegasus_vlm::solver::{FlightCondition, VortexLatticeSolver};

const NX: usize = 2;
const NY: usize = 3;

/// A swept, tapered, dihedralled mesh.
///
/// Deliberately not rectangular: a flat untwisted panel makes several geometry
/// terms vanish at once, so it can hide a sign error in the normal or the moment
/// arm that a general planform exposes.
fn mesh_input(duplicated: bool) -> MeshInput {
    let mut surface = SurfaceMesh::new("Wing", NX, NY);
    for j in 0..=NY {
        let fraction = j as f64 / NY as f64;
        let y = 2.0 * fraction;
        let leading_edge = 0.35 * y;
        let chord = 1.0 - 0.4 * fraction;
        let z = 0.15 * y;
        for i in 0..=NX {
            let x = leading_edge + chord * i as f64 / NX as f64;
            surface.set_node(i, j, Vec3::new(x, y, z)).unwrap();
        }
    }
    if duplicated {
        surface.y_duplicated = true;
        surface.y_duplicate_plane = 0.0;
    }
    let mut input = MeshInput::default();
    input.surfaces.push(surface);
    input.reference.area = 1.6;
    input.reference.chord = 0.8;
    input.reference.span = 2.0;
    input.reference.moment_reference_point = Vec3::new(0.25, 0.0, 0.0);
    input
}

fn condition() -> FlightCondition {
    FlightCondition {
        alpha: 4.0 * DEGREES_TO_RADIANS,
        beta: 2.0 * DEGREES_TO_RADIANS,
        roll_rate: 0.03,
        pitch_rate: -0.02,
        yaw_rate: 0.01,
        ..FlightCondition::default()
    }
}

/// Moment points that do not move with the mesh, matching how a coupled solver
/// supplies its beam control points.
fn moment_points(strip_count: usize) -> Vec<Vec3> {
    (0..strip_count)
        .map(|strip| {
            let fraction = (strip as f64 + 0.5) / strip_count as f64;
            Vec3::new(0.3 + 0.1 * fraction, 2.0 * fraction - 1.0, 0.05)
        })
        .collect()
}

/// Two independent functionals of the strip loads, with no symmetry that could
/// mask a transposed index.
fn cotangents(strip_count: usize) -> Vec<StripLoadCotangent> {
    (0..2)
        .map(|column| {
            let mut cotangent = StripLoadCotangent::zeros(strip_count);
            for strip in 0..strip_count {
                let seed = (strip + 1) as f64 + 0.5 * column as f64;
                cotangent.force[strip] = Vec3::new(0.7 * seed, -0.3 * seed, 1.1 / seed);
                cotangent.moment[strip] = Vec3::new(-0.4 / seed, 0.9 * seed, 0.25 * seed);
            }
            cotangent
        })
        .collect()
}

/// Solves the mesh and evaluates every functional at that solution.
///
/// The resultants are rebuilt from `StripLoads` exactly as a coupled solver
/// does: the reported force and moment are area- and chord-normalized
/// coefficients about the strip quarter chord, carrying AVL's roll and yaw sign
/// flips, so they are redimensionalized, unflipped, and transported to the
/// functional's moment point.
fn functionals(input: &MeshInput, points: &[Vec3], columns: &[StripLoadCotangent]) -> Vec<f64> {
    functionals_at(input, points, columns, &condition())
}

/// The same, at a flight condition the caller chooses.
fn functionals_at(
    input: &MeshInput,
    points: &[Vec3],
    columns: &[StripLoadCotangent],
    condition: &FlightCondition,
) -> Vec<f64> {
    let mut solver = VortexLatticeSolver::new();
    solver.set_mesh(input).unwrap();
    let solution = solver.solve(condition).unwrap();
    let results = ForceCalculator::new().compute_with_references(
        &input.reference,
        solver.mesh(),
        &solution,
        condition,
    );

    let mut values = vec![0.0; columns.len()];
    for surface in &results.surface_loads {
        for loads in &surface.strip_loads {
            let force = loads.body_force * loads.area;
            let moment = Vec3::new(
                -loads.body_moment.x,
                loads.body_moment.y,
                -loads.body_moment.z,
            ) * loads.chord
                * loads.area
                + (loads.quarter_chord_position - points[loads.strip_index]).cross(&force);
            for (value, column) in values.iter_mut().zip(columns) {
                *value += column.force[loads.strip_index].dot(&force)
                    + column.moment[loads.strip_index].dot(&moment);
            }
        }
    }
    values
}

/// Differences every node coordinate of the surface against the adjoint.
fn check(duplicated: bool) {
    let input = mesh_input(duplicated);
    let mut solver = VortexLatticeSolver::new();
    solver.set_mesh(&input).unwrap();
    solver.solve(&condition()).unwrap();

    let strip_count = solver.mesh().strip_count();
    assert_eq!(strip_count, if duplicated { 2 * NY } else { NY });
    let points = moment_points(strip_count);
    let columns = cotangents(strip_count);

    let condition = condition();
    let adjoint = StripLoadAdjoint::new(&solver, &condition, &points).unwrap();
    let cotangent = adjoint.mesh_node_cotangents(&input, &columns).unwrap();
    assert_eq!(cotangent.len(), columns.len());

    let step = 1.0e-7;
    let mut worst = 0.0f64;
    for node in 0..input.surfaces[0].node_count() {
        for direction in 0..3 {
            let mut plus = input.clone();
            let mut minus = input.clone();
            plus.surfaces[0].nodes[node][direction] += step;
            minus.surfaces[0].nodes[node][direction] -= step;
            let high = functionals(&plus, &points, &columns);
            let low = functionals(&minus, &points, &columns);

            for column in 0..columns.len() {
                let difference = (high[column] - low[column]) / (2.0 * step);
                let analytical = cotangent[column].surfaces[0][node][direction];
                let scale = 1.0 + difference.abs().max(analytical.abs());
                let error = (analytical - difference).abs() / scale;
                worst = worst.max(error);
                assert!(
                    error < 1.0e-6,
                    "column {column}, node {node}, direction {direction}: \
                     adjoint {analytical:e} vs difference {difference:e}"
                );
            }
        }
    }
    assert!(worst > 0.0, "every component was identically zero");
}

#[test]
fn mesh_node_cotangents_match_finite_differences() {
    check(false);
}

#[test]
fn mesh_node_cotangents_match_finite_differences_on_a_duplicated_surface() {
    check(true);
}

#[test]
fn a_symmetry_plane_is_rejected_rather_than_silently_ignored() {
    use pegasus_vlm::core::SymmetryType;
    use pegasus_vlm::solver::SymmetryConfig;

    let input = mesh_input(false);
    let mut solver = VortexLatticeSolver::new();
    solver.set_mesh(&input).unwrap();
    solver.set_symmetry(SymmetryConfig {
        y_symmetry: SymmetryType::Symmetric,
        ..SymmetryConfig::default()
    });
    solver.solve(&condition()).unwrap();

    let condition = condition();
    let points = moment_points(solver.mesh().strip_count());
    assert!(StripLoadAdjoint::new(&solver, &condition, &points).is_err());
}

#[test]
fn a_moment_point_count_mismatch_is_an_error() {
    let input = mesh_input(false);
    let mut solver = VortexLatticeSolver::new();
    solver.set_mesh(&input).unwrap();
    solver.solve(&condition()).unwrap();

    let condition = condition();
    assert!(StripLoadAdjoint::new(&solver, &condition, &[Vec3::zeros()]).is_err());
}

/// The flight-condition half of the same adjoint, against central differences.
///
/// The angles reach a strip's load explicitly, through the total velocity the
/// Kutta-Joukowski force is built from, and implicitly, through the
/// flow-tangency right-hand side that sets the circulation. A check on only one
/// of them would pass with the other missing, so this differences the whole
/// functional and both angles.
#[test]
fn flight_condition_cotangents_match_finite_differences() {
    let input = mesh_input(false);
    let mut solver = VortexLatticeSolver::new();
    solver.set_mesh(&input).unwrap();
    solver.solve(&condition()).unwrap();

    let strip_count = solver.mesh().strip_count();
    let points = moment_points(strip_count);
    let columns = cotangents(strip_count);

    let flight = condition();
    let adjoint = StripLoadAdjoint::new(&solver, &flight, &points).unwrap();
    let derivatives = adjoint.flight_condition_cotangents(&columns).unwrap();
    assert_eq!(derivatives.len(), columns.len());

    let step = 1e-7;
    for (name, is_alpha) in [("alpha", true), ("beta", false)] {
        let mut evaluated = Vec::with_capacity(2);
        for sign in [1.0, -1.0] {
            let mut perturbed = condition();
            if is_alpha {
                perturbed.alpha += sign * step;
            } else {
                perturbed.beta += sign * step;
            }
            evaluated.push(functionals_at(&input, &points, &columns, &perturbed));
        }

        for (column, total) in derivatives.iter().enumerate() {
            let fd = (evaluated[0][column] - evaluated[1][column]) / (2.0 * step);
            let analytical = if is_alpha {
                total.by_alpha
            } else {
                total.by_beta
            };
            let scale = 1.0_f64.max(analytical.abs()).max(fd.abs());
            assert!(
                (analytical - fd).abs() / scale < 1e-6,
                "d(functional {column})/d({name}): {analytical:e} vs {fd:e}"
            );
        }
    }
}

/// A cotangent that does not describe this lattice is an error, not a wrong
/// answer.
#[test]
fn a_flight_condition_cotangent_of_the_wrong_length_is_an_error() {
    let input = mesh_input(false);
    let mut solver = VortexLatticeSolver::new();
    solver.set_mesh(&input).unwrap();
    solver.solve(&condition()).unwrap();
    let points = moment_points(solver.mesh().strip_count());
    let flight = condition();
    let adjoint = StripLoadAdjoint::new(&solver, &flight, &points).unwrap();

    let short = vec![StripLoadCotangent::zeros(solver.mesh().strip_count() - 1)];
    assert!(adjoint.flight_condition_cotangents(&short).is_err());
    assert!(adjoint.flight_condition_cotangents(&[]).unwrap().is_empty());
}
