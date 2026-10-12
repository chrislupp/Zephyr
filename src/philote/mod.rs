//! Philote-MDO gRPC disciplines exposing the solver to an MDO framework.
//!
//! Enabled by the non-default `philote` feature. Two explicit disciplines are
//! provided, and they chain: [`MeshDiscipline`] turns a section-based
//! parameterization into a node mesh plus reference quantities, and
//! [`VlmDiscipline`] takes that mesh and returns the six force and moment
//! coefficients, optionally dimensional loads and solved trim variables.
//! Internal trim keeps the service explicit, with analytical derivatives through
//! the converged operating condition.
//!
//! Both supply analytical partials. The VLM discipline's node partials come
//! from the adjoint, so a gradient-based optimizer pays one transposed solve
//! per output rather than one primal solve per design variable.
//!
//! # Options
//!
//! Philote delivers options as JSON. Per-surface lists accept either a JSON
//! array — `"surface_nx": [8, 4]` — or the comma-separated string the v0.2.0
//! C++ used, `"surface_nx": "8,4"`. Both are accepted so a client written
//! against either implementation works unchanged.
//!
//! # Example
//!
//! ```no_run
//! use pegasus_vlm::philote::VlmDiscipline;
//! use philote_mdo::server::ExplicitServer;
//! # fn main() {
//! let server = ExplicitServer::new(VlmDiscipline::default());
//! # let _ = server;
//! # }
//! ```
//!
//! See `examples/philote_vlm_server.rs` and `examples/philote_mesh_server.rs`
//! for runnable servers.

use ndarray::{ArrayD, IxDyn};
use philote_mdo::{PhiloteError, Result};

mod aero_evaluation;
mod mesh;
mod trim;
mod vlm;

pub use mesh::MeshDiscipline;
pub use vlm::{SurfaceOption, VlmDiscipline};

/// A configuration or option problem.
pub(crate) fn invalid(message: impl Into<String>) -> PhiloteError {
    PhiloteError::ConfigurationError(message.into())
}

/// A solver refusal, tagged with the discipline call that produced it so a
/// client can tell a bad mesh from a bad flight condition.
pub(crate) fn rejected(context: &str, message: impl Into<String>) -> PhiloteError {
    PhiloteError::validation(context, message)
}

/// A one-dimensional array of the given length.
pub(crate) fn flat(values: Vec<f64>) -> ArrayD<f64> {
    let length = values.len();
    ArrayD::from_shape_vec(IxDyn(&[length]), values).expect("length matches the shape")
}

/// A single-element array, Philote's spelling of a scalar.
pub(crate) fn scalar(value: f64) -> ArrayD<f64> {
    flat(vec![value])
}

/// Read a required input, reporting the name if it is absent.
pub(crate) fn input<'a>(inputs: &'a philote_mdo::ArrayMap, name: &str) -> Result<&'a ArrayD<f64>> {
    inputs
        .get(name)
        .ok_or_else(|| PhiloteError::VariableNotFound(name.to_string()))
}

/// Read a required scalar input.
pub(crate) fn scalar_input(inputs: &philote_mdo::ArrayMap, name: &str) -> Result<f64> {
    let array = input(inputs, name)?;
    array
        .iter()
        .next()
        .copied()
        .ok_or_else(|| PhiloteError::ShapeMismatch {
            expected: vec![1],
            actual: array.shape().to_vec(),
        })
}

/// Read an input expected to hold exactly `expected` values.
pub(crate) fn sized_input(
    inputs: &philote_mdo::ArrayMap,
    name: &str,
    expected: usize,
) -> Result<Vec<f64>> {
    let array = input(inputs, name)?;
    if array.len() != expected {
        return Err(PhiloteError::ShapeMismatch {
            expected: vec![expected],
            actual: array.shape().to_vec(),
        });
    }
    Ok(array.iter().copied().collect())
}

/// Split a per-surface option that may arrive as a JSON array or as the
/// comma-separated string the C++ discipline used.
///
/// Returning the raw tokens rather than parsed numbers lets one helper serve
/// the integer, float and boolean lists alike.
pub(crate) fn option_tokens(value: &serde_json::Value, name: &str) -> Result<Vec<String>> {
    match value {
        serde_json::Value::Array(items) => items
            .iter()
            .map(|item| match item {
                serde_json::Value::String(text) => Ok(text.trim().to_string()),
                other => Ok(other.to_string()),
            })
            .collect(),
        serde_json::Value::String(text) => Ok(text
            .split(',')
            .map(|token| token.trim().to_string())
            .filter(|token| !token.is_empty())
            .collect()),
        other => Err(invalid(format!(
            "option `{name}` must be an array or a comma-separated string, got {other}"
        ))),
    }
}

/// Parse a token as a boolean, accepting the spellings AVL-style inputs use.
pub(crate) fn parse_bool(token: &str) -> bool {
    matches!(
        token.trim().to_ascii_lowercase().as_str(),
        "true" | "1" | "yes" | "y"
    )
}

/// Read a scalar option as an integer.
pub(crate) fn option_usize(value: &serde_json::Value, name: &str) -> Result<usize> {
    let parsed = match value {
        serde_json::Value::Number(number) => number.as_f64(),
        serde_json::Value::String(text) => text.trim().parse::<f64>().ok(),
        _ => None,
    };
    let parsed = parsed.ok_or_else(|| invalid(format!("option `{name}` must be a number")))?;
    if !parsed.is_finite() || parsed < 0.0 {
        return Err(invalid(format!(
            "option `{name}` must be a non-negative number, got {parsed}"
        )));
    }
    Ok(parsed as usize)
}

/// Read a scalar option as a boolean.
pub(crate) fn option_bool(value: &serde_json::Value) -> bool {
    match value {
        serde_json::Value::Bool(flag) => *flag,
        serde_json::Value::String(text) => parse_bool(text),
        serde_json::Value::Number(number) => number.as_f64().is_some_and(|value| value != 0.0),
        _ => false,
    }
}
