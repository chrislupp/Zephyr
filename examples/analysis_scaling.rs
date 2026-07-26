use std::time::Instant;

use pegasus_vlm::core::{DEGREES_TO_RADIANS, Vec3};
use pegasus_vlm::geometry::{MeshInput, ReferenceQuantities, SurfaceMesh};
use pegasus_vlm::solver::{FlightCondition, VortexLatticeSolver};

fn mesh(ny: usize) -> MeshInput {
    let nx = 6;
    let mut surface = SurfaceMesh::new("wing", nx, ny);
    surface.y_duplicated = true;
    for i in 0..=nx {
        let chord_fraction = 0.5 * (1.0 - (std::f64::consts::PI * i as f64 / nx as f64).cos());
        for j in 0..=ny {
            surface
                .set_node(
                    i,
                    j,
                    Vec3::new(chord_fraction, 5.0 * j as f64 / ny as f64, 0.0),
                )
                .expect("indices are in range");
        }
    }
    MeshInput {
        surfaces: vec![surface],
        reference: ReferenceQuantities {
            area: 10.0,
            chord: 1.0,
            span: 10.0,
            moment_reference_point: Vec3::new(0.25, 0.0, 0.0),
        },
    }
}

fn main() -> pegasus_vlm::Result<()> {
    let condition = FlightCondition {
        alpha: 5.0 * DEGREES_TO_RADIANS,
        dynamic_pressure: 0.5 * 1.225 * 25.0_f64.powi(2),
        ..FlightCondition::default()
    };
    println!("ny,panels,time_ms");
    for ny in [8, 12, 16, 24, 32, 48, 64, 96] {
        let input = mesh(ny);
        let repetitions = 10;
        let start = Instant::now();
        for _ in 0..repetitions {
            let mut solver = VortexLatticeSolver::new();
            solver.set_mesh(&input)?;
            solver.solve(&condition)?;
        }
        let milliseconds = start.elapsed().as_secs_f64() * 1_000.0 / repetitions as f64;
        println!("{ny},{},{milliseconds:.6}", 6 * 2 * ny);
    }
    Ok(())
}
