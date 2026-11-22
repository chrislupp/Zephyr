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

//! Error types returned by Zephyr.

use std::path::PathBuf;

/// Crate result type.
pub type Result<T> = std::result::Result<T, ZephyrError>;

/// Where a parse failure occurred.
///
/// Line and path are separate fields so a caller can act on them — jump to the
/// line, group errors by file — rather than parsing them back out of a
/// rendered message.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct ParseLocation {
    /// Source file, when the input came from one.
    pub path: Option<PathBuf>,
    /// 1-based line number, when known.
    pub line: Option<usize>,
}

impl std::fmt::Display for ParseLocation {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match (&self.path, self.line) {
            (Some(path), Some(line)) => write!(formatter, " in {} at line {line}", path.display()),
            (Some(path), None) => write!(formatter, " in {}", path.display()),
            (None, Some(line)) => write!(formatter, " at line {line}"),
            (None, None) => Ok(()),
        }
    }
}

/// Failures produced while constructing, solving, parsing, or exporting.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum ZephyrError {
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
        /// Structured path/line location.
        location: ParseLocation,
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
    /// Invalid trim configuration or prescribed trim data.
    #[error("invalid trim problem: {0}")]
    InvalidTrim(String),
    /// A trim solve failed to reach a locally differentiable equilibrium.
    #[error("trim failed: {0}")]
    TrimFailed(String),
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

impl ZephyrError {
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
            location: ParseLocation {
                path: None,
                line: Some(line),
            },
            message: message.into(),
        }
    }

    /// Construct a path-input parse error.
    #[must_use]
    pub fn parse_path(path: &std::path::Path, line: usize, message: impl Into<String>) -> Self {
        Self::Parse {
            location: ParseLocation {
                path: Some(path.to_owned()),
                line: Some(line),
            },
            message: message.into(),
        }
    }

    /// Attach a source path to a parse error that was produced without one.
    ///
    /// Used when text is parsed by a string entry point on behalf of a file
    /// entry point, so the file is named in the message without the inner
    /// parser needing to know it.
    #[must_use]
    pub fn with_source_path(self, path: impl Into<PathBuf>) -> Self {
        match self {
            Self::Parse {
                mut location,
                message,
            } => {
                location.path.get_or_insert_with(|| path.into());
                Self::Parse { location, message }
            }
            other => other,
        }
    }
}
