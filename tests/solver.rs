mod common;

use approx::assert_abs_diff_eq;
use pegasus_vlm::core::{DEGREES_TO_RADIANS, SpacingType, SymmetryType, Vec3, VecX};
use pegasus_vlm::geometry::{Aircraft, Section, Surface};
use pegasus_vlm::solver::{
    FlightCondition, InfluenceMatrixBuilder, SymmetryConfig, VortexLattice, VortexLatticeSolver,
    horseshoe_vortex_velocity,
};

use common::{rectangular_mesh, rectangular_wing};

#[test]
fn aircraft_lattice_preserves_expected_panel_order_and_counts() {
    let aircraft = rectangular_wing(3, 5);
    let mut lattice = VortexLattice::default();
    lattice.generate_from_aircraft(&aircraft).unwrap();
    assert_eq!(lattice.vortex_count(), 15);
    assert_eq!(lattice.strip_count(), 5);
    assert_eq!(lattice.surface_count(), 1);
    for (strip_index, strip) in lattice.strips().iter().enumerate() {
        assert_eq!(strip.first_vortex_index, strip_index * 3);
        assert_eq!(strip.num_vortices, 3);
        for offset in 0..3 {
            assert_eq!(
                lattice.vortices()[strip.first_vortex_index + offset].strip_index,
                strip_index
            );
        }
    }
}

#[test]
fn flat_rectangular_lattice_has_quarter_chord_vortices_and_three_quarter_controls() {
    let aircraft = rectangular_wing(1, 2);
    let mut lattice = VortexLattice::default();
    lattice.generate_from_aircraft(&aircraft).unwrap();
    for vortex in lattice.vortices() {
        assert_abs_diff_eq!(vortex.bound_vortex_center.x, 0.25, epsilon = 1.0e-12);
        assert_abs_diff_eq!(vortex.control_point.x, 0.75, epsilon = 1.0e-12);
        assert_abs_diff_eq!(vortex.normal_at_control_point, Vec3::z(), epsilon = 1.0e-12);
    }
}

#[test]
fn twist_rotates_panel_normals() {
    let mut wing = Surface::new("Twisted");
    wing.set_chordwise_panel_count(1);
    wing.set_spanwise_panel_count(2);
    wing.add_section(Section::new(Vec3::new(0.0, 0.0, 0.0), 1.0, 10.0));
    wing.add_section(Section::new(Vec3::new(0.0, 2.0, 0.0), 1.0, 10.0));
    let mut aircraft = Aircraft::new("Twisted");
    aircraft.add_surface(wing);
    let mut lattice = VortexLattice::default();
    lattice.generate_from_aircraft(&aircraft).unwrap();
    assert!(lattice.vortices()[0].normal_at_control_point.x > 0.1);
    assert_abs_diff_eq!(
        lattice.vortices()[0].normal_at_control_point.norm(),
        1.0,
        epsilon = 1.0e-12
    );
}

#[test]
fn y_duplicate_mirrors_geometry_and_preserves_orientation() {
    let mut aircraft = rectangular_wing(2, 3);
    let wing = &mut aircraft.surfaces_mut()[0];
    wing.set_y_duplicate(0.0);
    let mut lattice = VortexLattice::default();
    lattice.generate_from_aircraft(&aircraft).unwrap();
    assert_eq!(lattice.surface_count(), 2);
    assert_eq!(lattice.vortex_count(), 12);
    let original = &lattice.surface_info()[0];
    let mirror = &lattice.surface_info()[1];
    assert!(mirror.is_y_duplicated);
    for offset in 0..original.num_vortices {
        let left = &lattice.vortices()[original.first_vortex_index + offset];
        let right = &lattice.vortices()[mirror.first_vortex_index + offset];
        assert_abs_diff_eq!(left.control_point.x, right.control_point.x);
        assert_abs_diff_eq!(left.control_point.y, -right.control_point.y);
        assert!(right.bound_vortex_right.y > right.bound_vortex_left.y);
    }
}

#[test]
fn direct_mesh_and_aircraft_paths_match_for_rectangular_geometry() {
    let mut aircraft = rectangular_wing(2, 4);
    aircraft.surfaces_mut()[0].set_chordwise_spacing(SpacingType::Uniform);
    aircraft.surfaces_mut()[0].set_spanwise_spacing(SpacingType::Uniform);
    let mesh = rectangular_mesh(2, 4);
    let mut aircraft_solver = VortexLatticeSolver::new();
    let mut mesh_solver = VortexLatticeSolver::new();
    aircraft_solver.set_aircraft(&aircraft).unwrap();
    mesh_solver.set_mesh(&mesh).unwrap();
    let condition = FlightCondition {
        alpha: 4.0 * DEGREES_TO_RADIANS,
        ..FlightCondition::default()
    };
    let aircraft_solution = aircraft_solver.solve(&condition).unwrap();
    let mesh_solution = mesh_solver.solve(&condition).unwrap();
    assert_abs_diff_eq!(
        aircraft_solution.circulation,
        mesh_solution.circulation,
        epsilon = 1.0e-11
    );
}

#[test]
fn direct_mesh_dimensions_and_normal_recalculation_are_valid() {
    let mesh = rectangular_mesh(2, 3);
    let mut lattice = VortexLattice::default();
    lattice.generate_from_mesh(&mesh).unwrap();
    assert_eq!(lattice.vortex_count(), 6);
    assert_eq!(lattice.strip_count(), 3);
    assert!(lattice.mesh_input().is_some());
    lattice.recalculate_normals();
    assert!(
        lattice
            .vortices()
            .iter()
            .all(|vortex| (vortex.normal_at_control_point.norm() - 1.0).abs() < 1.0e-12)
    );
}

#[test]
fn regularized_horseshoe_kernel_is_finite_near_the_vortex() {
    let left = Vec3::new(0.0, -0.5, 0.0);
    let right = Vec3::new(0.0, 0.5, 0.0);
    for point in [
        Vec3::new(0.0, 0.0, 1.0e-10),
        Vec3::new(1.0, 0.0, 0.2),
        Vec3::new(-1.0, 0.2, 0.1),
    ] {
        let velocity = horseshoe_vortex_velocity(&point, &left, &right, 0.9, 0.01, true);
        assert!(velocity.iter().all(|value| value.is_finite()));
    }
}

#[test]
fn excluding_bound_segment_changes_only_the_finite_segment_contribution() {
    let left = Vec3::new(0.0, -0.5, 0.0);
    let right = Vec3::new(0.0, 0.5, 0.0);
    let point = Vec3::new(0.5, 0.0, 0.2);
    let with_bound = horseshoe_vortex_velocity(&point, &left, &right, 1.0, 0.001, true);
    let without_bound = horseshoe_vortex_velocity(&point, &left, &right, 1.0, 0.001, false);
    assert!((with_bound - without_bound).norm() > 1.0e-4);
}

#[test]
fn influence_matrix_dimensions_and_normal_projection_are_consistent() {
    let aircraft = rectangular_wing(2, 4);
    let mut lattice = VortexLattice::default();
    lattice.generate_from_aircraft(&aircraft).unwrap();
    let builder = InfluenceMatrixBuilder::new(&lattice);
    let (aic, velocity) = builder.normalwash_and_velocity();
    assert_eq!(aic.shape(), (8, 8));
    assert_eq!(velocity.shape(), (24, 8));
    for row in 0..8 {
        for column in 0..8 {
            let normal = lattice.vortices()[row].normal_at_control_point;
            let projected = normal.dot(&Vec3::new(
                velocity[(3 * row, column)],
                velocity[(3 * row + 1, column)],
                velocity[(3 * row + 2, column)],
            ));
            assert_abs_diff_eq!(aic[(row, column)], projected, epsilon = 1.0e-14);
        }
    }
}

#[test]
fn aic_is_nonsingular_and_unit_response_residual_is_small() {
    let aircraft = rectangular_wing(3, 6);
    let mut solver = VortexLatticeSolver::new();
    solver.set_aircraft(&aircraft).unwrap();
    solver.solve_unit_responses().unwrap();
    assert!(solver.is_aic_valid());
    assert!(solver.are_unit_responses_valid());
    let solution = solver.solve(&FlightCondition::default()).unwrap();
    let rhs = VecX::from_iterator(
        solver.mesh().vortex_count(),
        solver
            .mesh()
            .vortices()
            .iter()
            .map(|vortex| -vortex.normal_at_control_point.z),
    );
    let residual = solver.aic() * &solution.circulation_per_w - rhs;
    assert!(residual.norm() < 1.0e-10);
}

#[test]
fn zero_alpha_flat_wing_has_zero_circulation() {
    let aircraft = rectangular_wing(3, 6);
    let mut solver = VortexLatticeSolver::new();
    solver.set_aircraft(&aircraft).unwrap();
    let solution = solver.solve(&FlightCondition::default()).unwrap();
    assert!(solution.circulation.norm() < 1.0e-13);
    assert!(solution.induced_velocity.norm() < 1.0e-13);
}

#[test]
fn positive_and_negative_alpha_reverse_circulation() {
    let aircraft = rectangular_wing(3, 6);
    let solve = |alpha| {
        let mut solver = VortexLatticeSolver::new();
        solver.set_aircraft(&aircraft).unwrap();
        solver
            .solve(&FlightCondition {
                alpha,
                ..FlightCondition::default()
            })
            .unwrap()
            .circulation
    };
    let positive = solve(3.0 * DEGREES_TO_RADIANS);
    let negative = solve(-3.0 * DEGREES_TO_RADIANS);
    assert_abs_diff_eq!(positive, -negative, epsilon = 1.0e-12);
}

#[test]
fn freestream_direction_handles_alpha_and_beta() {
    let direction = VortexLatticeSolver::freestream_direction(&FlightCondition {
        alpha: 10.0 * DEGREES_TO_RADIANS,
        beta: 5.0 * DEGREES_TO_RADIANS,
        ..FlightCondition::default()
    });
    assert_abs_diff_eq!(direction.norm(), 1.0, epsilon = 1.0e-14);
    assert!(direction.x > 0.9);
    assert!(direction.y > 0.0);
    assert!(direction.z > 0.0);
}

#[test]
fn mach_change_rebuilds_the_cached_aic() {
    let aircraft = rectangular_wing(2, 4);
    let mut solver = VortexLatticeSolver::new();
    solver.set_aircraft(&aircraft).unwrap();
    solver.solve(&FlightCondition::default()).unwrap();
    let incompressible = solver.aic().clone();
    solver
        .solve(&FlightCondition {
            mach: 0.5,
            ..FlightCondition::default()
        })
        .unwrap();
    assert_eq!(solver.current_mach(), 0.5);
    assert!((&incompressible - solver.aic()).norm() > 1.0e-5);
}

#[test]
fn supersonic_mach_uses_the_documented_clamp() {
    let aircraft = rectangular_wing(2, 4);
    let mut solver = VortexLatticeSolver::new();
    solver.set_aircraft(&aircraft).unwrap();
    solver
        .solve(&FlightCondition {
            mach: 1.5,
            ..FlightCondition::default()
        })
        .unwrap();
    assert!(solver.aic().iter().all(|value| value.is_finite()));
}

#[test]
fn symmetry_images_change_the_influence_matrix() {
    let aircraft = rectangular_wing(2, 4);
    let mut solver = VortexLatticeSolver::new();
    solver.set_aircraft(&aircraft).unwrap();
    solver.rebuild_aic().unwrap();
    let base = solver.aic().clone();
    solver.set_symmetry(SymmetryConfig {
        y_symmetry: SymmetryType::Symmetric,
        y_symmetry_plane: 0.0,
        ..SymmetryConfig::default()
    });
    solver.rebuild_aic().unwrap();
    assert!((&base - solver.aic()).norm() > 1.0e-6);
}

#[test]
fn mesh_mutation_invalidates_dependent_caches() {
    let aircraft = rectangular_wing(2, 4);
    let mut solver = VortexLatticeSolver::new();
    solver.set_aircraft(&aircraft).unwrap();
    solver.solve(&FlightCondition::default()).unwrap();
    assert!(solver.is_aic_valid());
    solver.mutate_mesh(|mesh| {
        mesh.set_vortex_core_radius_chord_fraction(0.02);
    });
    assert!(!solver.is_aic_valid());
    assert!(!solver.are_unit_responses_valid());
    assert!(!solver.solution().is_valid);
}

#[test]
fn invalid_flight_conditions_return_errors() {
    let aircraft = rectangular_wing(1, 2);
    let mut solver = VortexLatticeSolver::new();
    solver.set_aircraft(&aircraft).unwrap();
    for condition in [
        FlightCondition {
            alpha: f64::NAN,
            ..FlightCondition::default()
        },
        FlightCondition {
            mach: -0.1,
            ..FlightCondition::default()
        },
        FlightCondition {
            dynamic_pressure: -1.0,
            ..FlightCondition::default()
        },
        FlightCondition {
            control_deflections: vec![("bad".into(), f64::INFINITY)],
            ..FlightCondition::default()
        },
    ] {
        assert!(solver.solve(&condition).is_err());
    }
}

#[test]
fn empty_solver_rejects_solve_and_aic_build() {
    let mut solver = VortexLatticeSolver::new();
    assert!(solver.rebuild_aic().is_err());
    assert!(solver.solve(&FlightCondition::default()).is_err());
}
