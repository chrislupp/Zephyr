//! Strict JSON configuration for the internal explicit-discipline trim solve.

use super::invalid;
use crate::trim::{AeroQuantity, ControlMixer, MeshTrimProblem, TrimBalance, TrimVariable};
use philote_mdo::Result;
use serde_json::{Map, Value};
use std::collections::{BTreeMap, HashMap};

fn object<'a>(v: &'a Value, allowed: &[&str], context: &str) -> Result<&'a Map<String, Value>> {
    let o = v
        .as_object()
        .ok_or_else(|| invalid(format!("{context} must be an object")))?;
    for key in o.keys() {
        if !allowed.contains(&key.as_str()) {
            return Err(invalid(format!("unknown {context} field '{key}'")));
        }
    }
    Ok(o)
}

fn number(v: &Value) -> Result<f64> {
    v.as_f64()
        .filter(|n| n.is_finite())
        .ok_or_else(|| invalid("trim numbers must be finite"))
}

fn limits(v: &Value) -> Result<[f64; 2]> {
    let a = v
        .as_array()
        .filter(|a| a.len() == 2)
        .ok_or_else(|| invalid("trim bounds require two finite numbers"))?;
    Ok([number(&a[0])?, number(&a[1])?])
}

fn strings(v: &Value) -> Result<Vec<String>> {
    v.as_array()
        .ok_or_else(|| invalid("trim names must be an array"))?
        .iter()
        .map(|v| {
            v.as_str()
                .map(String::from)
                .ok_or_else(|| invalid("trim names must be strings"))
        })
        .collect()
}

/// Decode the proposed `trim` option, accepting an object or JSON string.
pub(super) fn parse(value: &Value) -> Result<Option<MeshTrimProblem>> {
    let decoded = if let Some(s) = value.as_str() {
        serde_json::from_str(s).map_err(|e| invalid(format!("invalid trim JSON: {e}")))?
    } else {
        value.clone()
    };
    let o = object(
        &decoded,
        &[
            "enabled",
            "variables",
            "balances",
            "initial_guess",
            "bounds",
            "variable_scales",
            "max_iterations",
            "absolute_tolerance",
            "relative_tolerance",
            "rank_tolerance",
            "max_backtracks",
            "mixer",
        ],
        "trim",
    )?;
    if let Some(v) = o.get("enabled")
        && !v
            .as_bool()
            .ok_or_else(|| invalid("trim enabled must be boolean"))?
    {
        return Ok(None);
    }
    let names = strings(
        o.get("variables")
            .ok_or_else(|| invalid("trim requires variables"))?,
    )?;
    let mut problem = MeshTrimProblem::default();
    for name in names {
        let guess = o
            .get("initial_guess")
            .and_then(|v| v.get(&name))
            .map(number)
            .transpose()?
            .unwrap_or(0.0);
        let mut variable = TrimVariable::new(&name, guess);
        if let Some(v) = o.get("bounds").and_then(|v| v.get(&name)) {
            variable.bounds = limits(v)?;
        }
        if let Some(v) = o.get("variable_scales").and_then(|v| v.get(&name)) {
            variable.scale = number(v)?;
        }
        problem.variables.push(variable);
    }
    for key in ["initial_guess", "bounds", "variable_scales"] {
        if let Some(v) = o.get(key) {
            let map = v
                .as_object()
                .ok_or_else(|| invalid(format!("trim {key} must be an object")))?;
            if map
                .keys()
                .any(|name| !problem.variables.iter().any(|v| v.name == *name))
            {
                return Err(invalid(format!("trim {key} contains an unknown variable")));
            }
        }
    }
    let balances = o
        .get("balances")
        .and_then(Value::as_array)
        .ok_or_else(|| invalid("trim requires a balances array"))?;
    for v in balances {
        let row = object(v, &["quantity", "target_input", "scale"], "balance")?;
        let quantity = row
            .get("quantity")
            .and_then(Value::as_str)
            .and_then(AeroQuantity::from_name)
            .ok_or_else(|| invalid("unknown trim balance quantity"))?;
        let mut balance = TrimBalance::new(
            quantity,
            row.get("scale").map(number).transpose()?.unwrap_or(1.0),
        );
        if let Some(v) = row.get("target_input") {
            balance.target_input = v
                .as_str()
                .ok_or_else(|| invalid("target_input must be a string"))?
                .into();
        }
        problem.balances.push(balance);
    }
    for (key, target) in [
        (
            "absolute_tolerance",
            &mut problem.options.absolute_tolerance,
        ),
        (
            "relative_tolerance",
            &mut problem.options.relative_tolerance,
        ),
        ("rank_tolerance", &mut problem.options.rank_tolerance),
    ] {
        if let Some(v) = o.get(key) {
            *target = number(v)?;
        }
    }
    for (key, target) in [
        ("max_iterations", &mut problem.options.max_iterations),
        ("max_backtracks", &mut problem.options.max_backtracks),
    ] {
        if let Some(v) = o.get(key) {
            *target = v
                .as_u64()
                .and_then(|n| usize::try_from(n).ok())
                .ok_or_else(|| invalid(format!("{key} must be a nonnegative integer")))?;
        }
    }
    if let Some(v) = o.get("mixer") {
        let m = object(
            v,
            &["commands", "surfaces", "matrix", "offsets", "bounds"],
            "mixer",
        )?;
        let get = |key| {
            m.get(key)
                .ok_or_else(|| invalid(format!("mixer requires {key}")))
        };
        let commands = strings(get("commands")?)?;
        let surfaces = strings(get("surfaces")?)?;
        let matrix = get("matrix")?
            .as_array()
            .ok_or_else(|| invalid("mixer matrix must be an array"))?
            .iter()
            .map(|v| {
                v.as_array()
                    .ok_or_else(|| invalid("mixer rows must be arrays"))?
                    .iter()
                    .map(number)
                    .collect()
            })
            .collect::<Result<_>>()?;
        let offsets = match m.get("offsets") {
            Some(v) => v
                .as_array()
                .ok_or_else(|| invalid("mixer offsets must be an array"))?
                .iter()
                .map(number)
                .collect::<Result<_>>()?,
            None => vec![0.0; surfaces.len()],
        };
        let mut bounds = BTreeMap::new();
        if let Some(v) = m.get("bounds") {
            for (name, limit) in v
                .as_object()
                .ok_or_else(|| invalid("mixer bounds must be an object"))?
            {
                bounds.insert(name.clone(), limits(limit)?);
            }
        }
        problem.mixer = Some(ControlMixer {
            commands,
            surfaces,
            matrix,
            offsets,
            bounds,
        });
    }
    Ok(Some(problem))
}

pub(super) fn targets(
    inputs: &philote_mdo::ArrayMap,
    problem: &MeshTrimProblem,
) -> Result<Vec<f64>> {
    problem
        .balances
        .iter()
        .map(|b| super::scalar_input(inputs, &b.target_input))
        .collect()
}

pub(super) fn external(
    inputs: &philote_mdo::ArrayMap,
    problem: &MeshTrimProblem,
) -> Result<[f64; 6]> {
    let mut loads = [0.0; 6];
    if problem.balances.iter().any(|b| b.quantity as usize >= 6) {
        for (i, q) in AeroQuantity::ALL[6..].iter().enumerate() {
            loads[i] = super::scalar_input(inputs, &format!("external_{}", q.name()))?;
        }
    }
    Ok(loads)
}

pub(super) fn selected(options: &HashMap<String, Value>) -> Result<Option<MeshTrimProblem>> {
    options
        .get("trim")
        .map(parse)
        .transpose()
        .map(Option::flatten)
}
