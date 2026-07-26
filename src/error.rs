//! Error types returned by Pegasus.

use std::path::PathBuf;

/// Crate result type.
pub type Result<T> = std::result::Result<T, PegasusError>;

/// Failures produced while constructing, solving, parsing, or exporting.
#[derive(Debug, thiserror::Error)]
pub enum PegasusError {
    /// Underlying filesystem error.
    #[error("I/O error for {path}: {source}")]
    Io {
        /// Affected path.
        path: PathBuf,
        /// Original error.
        #[source]
        source: std::io::Error,
    },
    /// Invalid text input.
    #[error("parse error{location}: {message}")]
    Parse {
        /// Human-readable path/line location.
        location: String,
        /// Error detail.
        message: String,
    },
    /// Invalid aircraft or airfoil definition.
    #[error("invalid geometry: {0}")]
    InvalidGeometry(String),
    /// Invalid structured mesh definition.
    #[error("invalid mesh: {0}")]
    InvalidMesh(String),
    /// Invalid flight condition.
    #[error("invalid flight condition: {0}")]
    InvalidFlightCondition(String),
    /// A dense solve could not produce a finite solution.
    #[error("singular or ill-conditioned aerodynamic system")]
    SingularSystem,
    /// Dynamic-mode decomposition failed.
    #[error("eigensolver failed: {0}")]
    Eigensolver(String),
    /// An operation requires a valid solved state.
    #[error("a valid aerodynamic solution is required")]
    InvalidSolution,
    /// A compile-time feature was disabled.
    #[error("feature '{0}' is disabled")]
    FeatureDisabled(&'static str),
}

impl PegasusError {
    /// Construct an I/O error carrying its path.
    #[must_use]
    pub fn io(path: impl Into<PathBuf>, source: std::io::Error) -> Self {
        Self::Io {
            path: path.into(),
            source,
        }
    }

    /// Construct a string-input parse error.
    #[must_use]
    pub fn parse(line: usize, message: impl Into<String>) -> Self {
        Self::Parse {
            location: format!(" at line {line}"),
            message: message.into(),
        }
    }

    /// Construct a path-input parse error.
    #[must_use]
    pub fn parse_path(path: &std::path::Path, line: usize, message: impl Into<String>) -> Self {
        Self::Parse {
            location: format!(" in {} at line {line}", path.display()),
            message: message.into(),
        }
    }
}
