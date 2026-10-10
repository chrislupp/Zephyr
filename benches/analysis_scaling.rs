use criterion::{BenchmarkId, Criterion, criterion_group, criterion_main};
use pegasus_vlm::core::{DEGREES_TO_RADIANS, Vec3};
use pegasus_vlm::geometry::{MeshInput, ReferenceQuantities, SurfaceMesh};
use pegasus_vlm::solver::{FlightCondition, VortexLatticeSolver};

fn mesh(ny: usize) -> MeshInput {
    let nx = 6;
    let mut surface = SurfaceMesh::new("wing", nx, ny);
    surface.y_duplicated = true;
    for i in 0..=nx {
        let x = 0.5 * (1.0 - (std::f64::consts::PI * i as f64 / nx as f64).cos());
        for j in 0..=ny {
            surface
                .set_node(i, j, Vec3::new(x, 5.0 * j as f64 / ny as f64, 0.0))
                .expect("indices are in range");
        }
    }
    MeshInput {
        consistent_forces: false,
        surfaces: vec![surface],
        reference: ReferenceQuantities {
            area: 10.0,
            chord: 1.0,
            span: 10.0,
            moment_reference_point: Vec3::new(0.25, 0.0, 0.0),
        },
    }
}

fn analysis_scaling(criterion: &mut Criterion) {
    let condition = FlightCondition {
        alpha: 5.0 * DEGREES_TO_RADIANS,
        ..FlightCondition::default()
    };
    let mut group = criterion.benchmark_group("analysis_scaling");
    for ny in [8, 12, 16, 24, 32, 48, 64, 96] {
        let input = mesh(ny);
        group.bench_with_input(
            BenchmarkId::new("panels", 12 * ny),
            &input,
            |bencher, mesh| {
                bencher.iter(|| {
                    let mut solver = VortexLatticeSolver::new();
                    solver.set_mesh(mesh).unwrap();
                    solver.solve(&condition).unwrap()
                });
            },
        );
    }
    group.finish();
}

criterion_group!(benches, analysis_scaling);
criterion_main!(benches);
