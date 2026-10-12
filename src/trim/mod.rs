//! Flight constraints and Newton trim.

mod constraints;
mod evaluation;
mod mesh_solver;
mod problem;
mod sensitivity;
mod solver;

pub use constraints::{
    ConstraintType, OperatingVariables, TrimConstraints, TrimPresets, VariableConstraint,
};
pub use evaluation::{AeroLinearization, AeroQuantity, aerodynamic_values, linearize_aerodynamics};
pub use mesh_solver::{MeshTrimSolution, MeshTrimSolver, TrimDiagnostics};
pub use problem::{ControlMixer, MeshTrimOptions, MeshTrimProblem, TrimBalance, TrimVariable};
pub use sensitivity::TrimPartials;
pub use solver::{TrimOptions, TrimResult, TrimSolver, TrimStatus};
