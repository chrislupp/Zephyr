//! Transpose solves and total flight-condition derivatives.

use crate::core::{MatX, Scalar, Vec3, VecX};
use crate::error::{PegasusError, Result};
use crate::solver::{FlightCondition, VortexLatticeSolver};

use super::{
    AeroOutput, force_alpha_gradient, force_beta_gradient, force_circulation_gradient,
    trefftz_circulation_gradients,
};

/// Adjoint solve backed by the primal AIC matrix.
pub struct AdjointSolver<'a> {
    solver: &'a VortexLatticeSolver,
}

impl<'a> AdjointSolver<'a> {
    pub fn new(solver: &'a VortexLatticeSolver) -> Result<Self> {
        if !solver.is_aic_valid() {
            return Err(PegasusError::InvalidSolution);
        }
        Ok(Self { solver })
    }

    pub fn solve(&self, rhs: &VecX) -> Result<VecX> {
        if rhs.len() != self.solver.aic().nrows() {
            return Err(PegasusError::InvalidSolution);
        }
        self.solver.solve_linearized_transpose(rhs)
    }

    pub fn solve_batched(&self, rhs: &MatX) -> Result<MatX> {
        if rhs.nrows() != self.solver.aic().nrows() {
            return Err(PegasusError::InvalidSolution);
        }
        self.solver.solve_linearized_transpose_batched(rhs)
    }

    #[must_use]
    pub fn unknown_count(&self) -> usize {
        self.solver.aic().nrows()
    }
}

/// Total derivatives for one selected coefficient.
#[derive(Clone, Copy, Debug, Default)]
pub struct TotalDerivative {
    pub by_alpha: Scalar,
    pub by_beta: Scalar,
}

impl TotalDerivative {
    /// Compute explicit plus implicit derivatives.
    pub fn compute(
        solver: &VortexLatticeSolver,
        condition: &FlightCondition,
        output: AeroOutput,
    ) -> Result<Self> {
        if !solver.solution().is_valid {
            return Err(PegasusError::InvalidSolution);
        }
        let gradient = if output == AeroOutput::InducedDrag {
            trefftz_circulation_gradients(
                solver.mesh(),
                solver.solution(),
                solver.reference_quantities().area,
                0.01,
                0.1,
            )
            .induced_drag
        } else {
            force_circulation_gradient(solver, condition, output)
        };
        let lambda = AdjointSolver::new(solver)?.solve(&gradient)?;
        let (sin_alpha, cos_alpha) = condition.alpha.sin_cos();
        let (sin_beta, cos_beta) = condition.beta.sin_cos();
        let velocity_alpha = Vec3::new(-sin_alpha * cos_beta, 0.0, cos_alpha * cos_beta);
        let velocity_beta = Vec3::new(-cos_alpha * sin_beta, cos_beta, -sin_alpha * sin_beta);
        let mut rhs_alpha = VecX::zeros(solver.mesh().vortex_count());
        let mut rhs_beta = VecX::zeros(solver.mesh().vortex_count());
        for (index, (vortex, normal)) in solver
            .mesh()
            .vortices()
            .iter()
            .zip(solver.linearized_normals())
            .enumerate()
        {
            if vortex.affected_by_freestream {
                rhs_alpha[index] = -normal.dot(&velocity_alpha);
                rhs_beta[index] = -normal.dot(&velocity_beta);
            }
        }
        let output_index = match output {
            AeroOutput::Lift => 0,
            AeroOutput::Drag => 1,
            AeroOutput::Side => 2,
            AeroOutput::Roll => 3,
            AeroOutput::Pitch => 4,
            AeroOutput::Yaw => 5,
            AeroOutput::InducedDrag => 0,
        };
        // These are the fixed-circulation explicit terms. In particular, the
        // alpha expression includes both freestream dependence and the
        // body-to-stability-axis rotation that was missing in v0.2.0.
        let (explicit_alpha, explicit_beta) = if output == AeroOutput::InducedDrag {
            (0.0, 0.0)
        } else {
            (
                force_alpha_gradient(
                    solver.mesh(),
                    solver.solution(),
                    condition,
                    solver.reference_quantities(),
                )[output_index],
                force_beta_gradient(
                    solver.mesh(),
                    solver.solution(),
                    condition,
                    solver.reference_quantities(),
                )[output_index],
            )
        };
        Ok(Self {
            by_alpha: explicit_alpha + lambda.dot(&rhs_alpha),
            by_beta: explicit_beta + lambda.dot(&rhs_beta),
        })
    }
}
