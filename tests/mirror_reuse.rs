/*
    Project Zephyr

    Copyright 2024-2026, Christopher A. Lupp

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This work has been cleared for public release, distribution unlimited, case
    number: AFRL-2026-2198.

    The views expressed are those of the authors and do not reflect the
    official guidance or position of the United States Government, the
    Department of Defense or of the United States Air Force.

    Statement from DoD: The Appearance of external hyperlinks does not
    constitute endorsement by the United States Department of Defense (DoD) of
    the linked websites, of the information, products, or services contained
    therein. The DoD does not exercise any editorial, security, or other
    control over the information you may find at these locations.
*/

//! Reflection reuse must be an optimization only: enabling it may never change
//! a single bit of output, and it must decline meshes it cannot represent
//! exactly rather than approximating them.

use zephyr_vlm::core::{Scalar, SymmetryType, Vec3};
use zephyr_vlm::geometry::{MeshInput, ReferenceQuantities, SurfaceMesh};
use zephyr_vlm::solver::{
    FlightCondition, InfluenceMatrixBuilder, SymmetryConfig, VortexLatticeSolver,
};

fn reference() -> ReferenceQuantities {
    ReferenceQuantities {
        area: 10.0,
        chord: 1.0,
        span: 10.0,
        moment_reference_point: Vec3::new(0.25, 0.0, 0.0),
    }
}

/// Flat y-duplicated wing, cosine-spaced chordwise.
fn duplicated_wing(nx: usize, ny: usize) -> MeshInput {
    let mut surface = SurfaceMesh::new("wing", nx, ny);
    surface.y_duplicated = true;
    for i in 0..=nx {
        let chord = 0.5 * (1.0 - (std::f64::consts::PI * i as Scalar / nx as Scalar).cos());
        for j in 0..=ny {
            surface
                .set_node(
                    i,
                    j,
                    Vec3::new(chord, 5.0 * j as Scalar / ny as Scalar, 0.0),
                )
                .unwrap();
        }
    }
    MeshInput {
        surfaces: vec![surface],
        reference: reference(),
    }
}

/// Swept and dihedralled, so panel normals have all three components non-zero
/// and the normalwash contraction is a genuine three-term sum.
fn oblique_duplicated(nx: usize, ny: usize) -> MeshInput {
    let mut surface = SurfaceMesh::new("swept", nx, ny);
    surface.y_duplicated = true;
    for i in 0..=nx {
        for j in 0..=ny {
            let span = j as Scalar / ny as Scalar;
            let chordwise = i as Scalar / nx as Scalar;
            surface
                .set_node(
                    i,
                    j,
                    Vec3::new(
                        0.6 * span + chordwise * (1.0 - 0.3 * span),
                        0.1 + 1.9 * span,
                        0.35 * span + chordwise * 0.15 * span,
                    ),
                )
                .unwrap();
        }
    }
    MeshInput {
        surfaces: vec![surface],
        reference: reference(),
    }
}

/// Duplicated about a plane that is not `y = 0`, so mirroring a coordinate is
/// `2*plane - y` and rounds. Reflection reuse must decline this.
fn offset_plane_wing(nx: usize, ny: usize) -> MeshInput {
    let mut surface = SurfaceMesh::new("wing", nx, ny);
    surface.y_duplicated = true;
    surface.y_duplicate_plane = 0.35;
    for i in 0..=nx {
        for j in 0..=ny {
            surface
                .set_node(
                    i,
                    j,
                    Vec3::new(
                        i as Scalar / nx as Scalar,
                        0.35 + 2.0 * j as Scalar / ny as Scalar,
                        0.0,
                    ),
                )
                .unwrap();
        }
    }
    MeshInput {
        surfaces: vec![surface],
        reference: reference(),
    }
}

/// Single surface with no duplicate, so no pairing exists.
fn plain_wing(nx: usize, ny: usize) -> MeshInput {
    let mut surface = SurfaceMesh::new("wing", nx, ny);
    for i in 0..=nx {
        for j in 0..=ny {
            surface
                .set_node(
                    i,
                    j,
                    Vec3::new(
                        i as Scalar / nx as Scalar,
                        -2.0 + 4.0 * j as Scalar / ny as Scalar,
                        0.0,
                    ),
                )
                .unwrap();
        }
    }
    MeshInput {
        surfaces: vec![surface],
        reference: reference(),
    }
}

fn condition(mach: Scalar) -> FlightCondition {
    FlightCondition {
        alpha: 0.06,
        beta: 0.03,
        roll_rate: 0.01,
        pitch_rate: 0.02,
        yaw_rate: 0.015,
        mach,
        dynamic_pressure: 0.5 * 1.225 * 25.0 * 25.0,
        ..FlightCondition::default()
    }
}

/// Every `f64` the solver produced, as raw bits, so comparison is exact rather
/// than tolerant.
fn output_bits(
    input: &MeshInput,
    symmetry: SymmetryConfig,
    mach: Scalar,
    mirror: bool,
) -> Vec<u64> {
    let mut solver = VortexLatticeSolver::new();
    solver.set_symmetry(symmetry);
    solver.set_mirror_reuse(mirror);
    solver.set_mesh(input).unwrap();
    let solution = solver.solve(&condition(mach)).unwrap();

    let mut bits = Vec::new();
    bits.extend(solver.aic().iter().map(|value| value.to_bits()));
    bits.extend(solver.velocity_matrix().iter().map(|value| value.to_bits()));
    for column in [
        &solution.circulation,
        &solution.circulation_per_u,
        &solution.circulation_per_v,
        &solution.circulation_per_w,
        &solution.circulation_per_p,
        &solution.circulation_per_q,
        &solution.circulation_per_r,
    ] {
        bits.extend(column.iter().map(|value| value.to_bits()));
    }
    bits.extend(
        solution
            .induced_velocity
            .iter()
            .map(|value| value.to_bits()),
    );
    bits
}

fn applies(input: &MeshInput, symmetry: SymmetryConfig, mach: Scalar) -> bool {
    let mut solver = VortexLatticeSolver::new();
    solver.set_symmetry(symmetry);
    solver.set_mesh(input).unwrap();
    let mut builder = InfluenceMatrixBuilder::new(solver.mesh());
    builder.set_mach_number(mach);
    builder.set_symmetry(symmetry);
    builder.set_mirror_reuse(true);
    builder.mirror_reuse_applies()
}

fn assert_bit_identical(label: &str, input: &MeshInput, symmetry: SymmetryConfig, mach: Scalar) {
    let plain = output_bits(input, symmetry, mach, false);
    let mirrored = output_bits(input, symmetry, mach, true);
    assert_eq!(
        plain.len(),
        mirrored.len(),
        "{label}: output length changed"
    );
    let differing = plain
        .iter()
        .zip(&mirrored)
        .filter(|(left, right)| left != right)
        .count();
    assert_eq!(
        differing,
        0,
        "{label}: reflection reuse changed {differing} of {} output values",
        plain.len()
    );
}

#[test]
fn reflection_reuse_is_bit_identical_where_it_engages() {
    let none = SymmetryConfig::default();
    for (label, input, mach) in [
        ("flat even span", duplicated_wing(6, 8), 0.0),
        ("flat odd span", duplicated_wing(6, 7), 0.0),
        ("flat compressible", duplicated_wing(6, 8), 0.2),
        ("oblique normals", oblique_duplicated(4, 6), 0.0),
        ("oblique compressible", oblique_duplicated(4, 6), 0.3),
        ("single chordwise row", duplicated_wing(1, 5), 0.0),
    ] {
        assert!(
            applies(&input, none, mach),
            "{label}: expected reflection reuse to engage, but it declined"
        );
        assert_bit_identical(label, &input, none, mach);
    }
}

#[test]
fn reflection_reuse_declines_meshes_it_cannot_represent_exactly() {
    let none = SymmetryConfig::default();
    let z_image = SymmetryConfig {
        z_symmetry: SymmetryType::Symmetric,
        z_symmetry_plane: 0.4,
        ..SymmetryConfig::default()
    };

    // A mirror plane away from the origin makes reflection a rounding
    // subtraction rather than a sign flip.
    assert!(!applies(&offset_plane_wing(4, 6), none, 0.0));
    // Symmetry images add terms whose equivariance is not established.
    assert!(!applies(&oblique_duplicated(4, 6), z_image, 0.0));
    // Nothing to pair.
    assert!(!applies(&plain_wing(4, 6), none, 0.0));

    // Declining must still produce correct results, not merely refuse.
    assert_bit_identical("offset plane", &offset_plane_wing(4, 6), none, 0.0);
    assert_bit_identical("z image", &oblique_duplicated(4, 6), z_image, 0.0);
    assert_bit_identical("not duplicated", &plain_wing(4, 6), none, 0.0);
}

/// Exercises the meshes large enough to reach the multithreaded assembly.
///
/// `drive_columns` hands off to rayon at `n >= PARALLEL_MIN_PANELS` (128 total
/// panels) and `drive_mirrored_columns` at `count >= PARALLEL_MIN_PANELS` (128
/// per half, so 256 total). Every other mesh in the suite is far below both
/// thresholds, so without this test `cargo test --all-features` compiles the
/// parallel code but never runs a single line of it, and CI's
/// `RAYON_NUM_THREADS` matrix sweeps a path no test reaches.
///
/// `duplicated_wing(8, 16)` is 128 panels per half and 256 total, which clears
/// both gates. Under `--features parallel` this is the parallel mirrored
/// assembly checked against the parallel plain assembly; under the default
/// build it is the same check on the serial drivers.
#[test]
fn large_meshes_reach_the_parallel_assembly_drivers() {
    let none = SymmetryConfig::default();

    // 128 per half / 256 total: clears both the plain and the mirrored gate.
    let both_gates = duplicated_wing(8, 16);
    assert!(applies(&both_gates, none, 0.0));
    assert_bit_identical("256 panels, both gates", &both_gates, none, 0.0);
    assert_bit_identical("256 panels, compressible", &both_gates, none, 0.25);

    // Oblique normals at scale, so the normalwash contraction is a genuine
    // three-term sum on the parallel path rather than a degenerate one.
    let oblique = oblique_duplicated(8, 16);
    assert!(applies(&oblique, none, 0.0));
    assert_bit_identical("256 panels, oblique", &oblique, none, 0.0);

    // Above the plain gate but below the mirrored one, so the two drivers
    // disagree about whether to distribute and must still agree bit-for-bit.
    let plain_gate_only = duplicated_wing(8, 9);
    assert_bit_identical("144 panels, plain gate only", &plain_gate_only, none, 0.0);

    // A non-duplicated mesh over the plain gate: reuse declines, so this is
    // `drive_columns` alone at n = 160.
    assert_bit_identical("160 panels, no duplicate", &plain_wing(8, 20), none, 0.0);
}

#[test]
fn reflection_reuse_leaves_the_default_path_alone() {
    // The flag defaults to off, and a solver that never touches it must behave
    // exactly like one explicitly opted out.
    let input = duplicated_wing(6, 8);
    let mut default_solver = VortexLatticeSolver::new();
    default_solver.set_mesh(&input).unwrap();
    assert!(!default_solver.mirror_reuse());

    let solution = default_solver.solve(&condition(0.0)).unwrap();
    let explicit = output_bits(&input, SymmetryConfig::default(), 0.0, false);
    let observed: Vec<u64> = default_solver
        .aic()
        .iter()
        .map(|value| value.to_bits())
        .collect();
    assert_eq!(observed, explicit[..observed.len()].to_vec());
    assert!(solution.is_valid);
}

/// Raw bits of every matrix the builder assembles, with reuse on or off.
fn builder_bits(
    input: &MeshInput,
    symmetry: SymmetryConfig,
    mach: Scalar,
    mirror: bool,
) -> (Vec<u64>, Vec<u64>, Vec<u64>) {
    let mut solver = VortexLatticeSolver::new();
    solver.set_mesh(input).unwrap();
    let mut builder = InfluenceMatrixBuilder::new(solver.mesh());
    builder.set_mach_number(mach);
    builder.set_symmetry(symmetry);
    builder.set_mirror_reuse(mirror);
    let to_bits = |matrix: &zephyr_vlm::core::MatX| -> Vec<u64> {
        matrix.iter().map(|value| value.to_bits()).collect()
    };
    let bound_only = builder.velocity_influence_matrix();
    let (aic, bound) = builder.normalwash_and_velocity();
    (to_bits(&bound_only), to_bits(&bound), to_bits(&aic))
}

/// The bound-centre (force-velocity) pass excludes each target's own bound
/// leg, a rule the control-point pass does not have. Reflection reuse must
/// still reproduce every one of its entries bit for bit, including the
/// diagonal and the mirror-image diagonal.
#[test]
fn bound_centre_reuse_is_bit_identical_entry_by_entry() {
    let none = SymmetryConfig::default();
    let z_image = SymmetryConfig {
        z_symmetry: SymmetryType::Symmetric,
        z_symmetry_plane: 0.4,
        ..SymmetryConfig::default()
    };
    for (label, input, symmetry, mach, engages) in [
        ("flat even span", duplicated_wing(6, 8), none, 0.0, true),
        ("flat odd span", duplicated_wing(6, 7), none, 0.0, true),
        ("flat compressible", duplicated_wing(6, 8), none, 0.2, true),
        ("oblique normals", oblique_duplicated(4, 6), none, 0.0, true),
        (
            "oblique compressible",
            oblique_duplicated(4, 6),
            none,
            0.3,
            true,
        ),
        (
            "single chordwise row",
            duplicated_wing(1, 5),
            none,
            0.0,
            true,
        ),
        // Over both parallel gates (256 panels, 128 per half).
        ("256 panels", duplicated_wing(8, 16), none, 0.0, true),
        (
            "256 panels oblique",
            oblique_duplicated(8, 16),
            none,
            0.25,
            true,
        ),
        // Declined: the plan refuses, so both runs evaluate every entry.
        ("offset plane", offset_plane_wing(4, 6), none, 0.0, false),
        ("z image", oblique_duplicated(4, 6), z_image, 0.0, false),
        ("not duplicated", plain_wing(4, 6), none, 0.0, false),
        (
            "not duplicated, 160 panels",
            plain_wing(8, 20),
            none,
            0.0,
            false,
        ),
    ] {
        assert_eq!(applies(&input, symmetry, mach), engages, "{label}");
        let plain = builder_bits(&input, symmetry, mach, false);
        let mirrored = builder_bits(&input, symmetry, mach, true);
        for (name, left, right) in [
            ("velocity_influence_matrix", &plain.0, &mirrored.0),
            ("normalwash_and_velocity bound", &plain.1, &mirrored.1),
            ("normalwash_and_velocity aic", &plain.2, &mirrored.2),
        ] {
            assert_eq!(left.len(), right.len(), "{label}: {name} shape changed");
            let differing = left.iter().zip(right).filter(|(a, b)| a != b).count();
            assert_eq!(
                differing,
                0,
                "{label}: {name} changed {differing} of {} entries",
                left.len()
            );
        }
        // Both passes of one builder agree with each other too.
        assert_eq!(plain.0, plain.1, "{label}: bound passes disagree");
    }
}
