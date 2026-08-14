//! Flight constraints and Newton trim.

mod constraints;
mod solver;

pub use constraints::{
    ConstraintType, OperatingVariables, TrimConstraints, TrimPresets, VariableConstraint,
};
pub use solver::{TrimOptions, TrimResult, TrimSolver, TrimStatus};
