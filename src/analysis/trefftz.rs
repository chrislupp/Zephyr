//! Trefftz integration and its analytical reverse, including reflected sources.
use super::operator::reflect;
use super::{Output, Work};
use crate::core::{EPSILON, PI, Scalar, Vec3, VecX};
use crate::geometry::ReferenceQuantities;
use crate::gradients::geometry_adjoint::LatticeCotangents;
use crate::solver::VortexLatticeSolver;
const K: Scalar = 1.0 / (2.0 * PI);

pub(crate) struct Far {
    pub lift: Scalar,
    pub side: Scalar,
    pub drag: Scalar,
}
impl Far {
    pub fn value(&self, o: Output, r: &ReferenceQuantities) -> Scalar {
        match o {
            Output::FarFieldLift => self.lift,
            Output::FarFieldSide => self.side,
            Output::FarFieldInducedDrag => self.drag,
            Output::SpanEfficiency => {
                if self.drag.abs() > EPSILON {
                    (self.lift * self.lift + self.side * self.side)
                        / (PI * r.span * r.span / r.area * self.drag)
                } else {
                    0.0
                }
            }
            _ => unreachable!(),
        }
    }
    pub fn weights(&self, o: Output, r: &ReferenceQuantities) -> [Scalar; 3] {
        match o {
            Output::FarFieldLift => [1.0, 0.0, 0.0],
            Output::FarFieldSide => [0.0, 1.0, 0.0],
            Output::FarFieldInducedDrag => [0.0, 0.0, 1.0],
            Output::SpanEfficiency => {
                if self.drag.abs() > EPSILON {
                    let scale = PI * r.span * r.span / r.area;
                    [
                        2.0 * self.lift / (scale * self.drag),
                        2.0 * self.side / (scale * self.drag),
                        -(self.lift * self.lift + self.side * self.side)
                            / (scale * self.drag * self.drag),
                    ]
                } else {
                    [0.0; 3]
                }
            }
            _ => unreachable!(),
        }
    }
}

struct Strip {
    left: Vec3,
    right: Vec3,
    gamma: Scalar,
    panel: usize,
    surface: usize,
    chord: Scalar,
}
fn strips(s: &VortexLatticeSolver) -> Vec<Strip> {
    s.mesh()
        .strips()
        .iter()
        .map(|strip| {
            let panel = strip.first_vortex_index + strip.num_vortices - 1;
            let v = &s.mesh().vortices()[panel];
            Strip {
                left: v.bound_vortex_left,
                right: v.bound_vortex_right,
                gamma: s
                    .solution()
                    .circulation
                    .rows(strip.first_vortex_index, strip.num_vortices)
                    .sum(),
                panel,
                surface: strip.surface_index,
                chord: strip.chord_center,
            }
        })
        .collect()
}

struct Pair {
    dy1: Scalar,
    dy2: Scalar,
    dz1: Scalar,
    dz2: Scalar,
    q1: Scalar,
    q2: Scalar,
    r1: Scalar,
    r2: Scalar,
    core: Scalar,
    chord_core: bool,
    width: Scalar,
    left: Vec3,
    right: Vec3,
    vy: Scalar,
    vz: Scalar,
}
fn pair(target: &Strip, source: &Strip, image: bool) -> Option<Pair> {
    let left = if image {
        reflect(source.right)
    } else {
        source.left
    };
    let right = if image {
        reflect(source.left)
    } else {
        source.right
    };
    let width = (right.y - left.y).hypot(right.z - left.z);
    let same = !image && target.surface == source.surface;
    let chord_core = 0.01 * source.chord >= 0.1 * width;
    let core = if same {
        0.0
    } else {
        (0.01 * source.chord).max(0.1 * width)
    };
    let center = 0.5 * (target.left + target.right);
    let dy1 = center.y - left.y;
    let dy2 = center.y - right.y;
    let dz1 = center.z - left.z;
    let dz2 = center.z - right.z;
    let q1 = dy1 * dy1 + dz1 * dz1;
    let q2 = dy2 * dy2 + dz2 * dz2;
    let r1 = (q1 * q1 + core.powi(4)).sqrt();
    let r2 = (q2 * q2 + core.powi(4)).sqrt();
    if r1 <= EPSILON || r2 <= EPSILON {
        return None;
    }
    Some(Pair {
        dy1,
        dy2,
        dz1,
        dz2,
        q1,
        q2,
        r1,
        r2,
        core,
        chord_core,
        width,
        left,
        right,
        vy: K * (dz1 / r1 - dz2 / r2),
        vz: K * (-dy1 / r1 + dy2 / r2),
    })
}

pub(crate) fn evaluate(s: &VortexLatticeSolver, half: bool, drag: bool, work: &mut Work) -> Far {
    let strips = strips(s);
    let area = s.reference_quantities().area;
    let mut result = Far {
        lift: 0.0,
        side: 0.0,
        drag: 0.0,
    };
    let scale = if half { 2.0 } else { 1.0 };
    for target in &strips {
        let dy = target.right.y - target.left.y;
        let dz = target.right.z - target.left.z;
        if dy.hypot(dz) < EPSILON {
            continue;
        }
        result.lift += scale * 2.0 * target.gamma * dy / area;
        if !half {
            result.side -= 2.0 * target.gamma * dz / area;
        }
        if drag {
            let (mut vy, mut vz) = (0.0, 0.0);
            for source in &strips {
                for image in 0..if half { 2 } else { 1 } {
                    work.trefftz_interactions += 1;
                    if let Some(p) = pair(target, source, image == 1) {
                        vy += source.gamma * p.vy;
                        vz += source.gamma * p.vz;
                    }
                }
            }
            result.drag += scale * target.gamma * (dz * vy - dy * vz) / area;
        }
    }
    result
}

/// Partial derivative with respect to circulation, and (when requested) geometry.
/// Only the requested scalar's reverse seed is accumulated.
pub(crate) fn reverse(
    s: &VortexLatticeSolver,
    half: bool,
    o: Output,
    geometry: bool,
    work: &mut Work,
) -> (VecX, LatticeCotangents) {
    let drag = matches!(o, Output::FarFieldInducedDrag | Output::SpanEfficiency);
    let far = if o == Output::SpanEfficiency {
        evaluate(s, half, true, work)
    } else {
        Far {
            lift: 0.0,
            side: 0.0,
            drag: 0.0,
        }
    };
    let [wl, ws, wd] = far.weights(o, s.reference_quantities());
    let scale = if half { 2.0 } else { 1.0 };
    let (wl, ws, wd) = (wl * scale, if half { 0.0 } else { ws }, wd * scale);
    let strips = strips(s);
    let area = s.reference_quantities().area;
    let mut gamma = vec![0.0; strips.len()];
    let mut cot = LatticeCotangents::zeros(s.mesh());
    for (i, target) in strips.iter().enumerate() {
        let dy = target.right.y - target.left.y;
        let dz = target.right.z - target.left.z;
        if dy.hypot(dz) < EPSILON {
            continue;
        }
        gamma[i] += 2.0 * (wl * dy - ws * dz) / area;
        if geometry {
            let v = Vec3::new(
                0.0,
                2.0 * wl * target.gamma / area,
                -2.0 * ws * target.gamma / area,
            );
            cot.bound_right[target.panel] += v;
            cot.bound_left[target.panel] -= v;
        }
        if !drag || wd == 0.0 {
            continue;
        }
        for (j, source) in strips.iter().enumerate() {
            for image in 0..if half { 2 } else { 1 } {
                work.trefftz_interactions += 1;
                let Some(p) = pair(target, source, image == 1) else {
                    continue;
                };
                let kernel = dz * p.vy - dy * p.vz;
                gamma[i] += wd * source.gamma * kernel / area;
                gamma[j] += wd * target.gamma * kernel / area;
                if !geometry {
                    continue;
                }
                let weight = wd * target.gamma * source.gamma / area;
                let by_width = Vec3::new(0.0, -weight * p.vz, weight * p.vy);
                cot.bound_right[target.panel] += by_width;
                cot.bound_left[target.panel] -= by_width;
                let sy = weight * dz;
                let sz = -weight * dy;
                let mut cy1 = -K * sz / p.r1;
                let mut cy2 = K * sz / p.r2;
                let mut cz1 = K * sy / p.r1;
                let mut cz2 = -K * sy / p.r2;
                let cr1 = K * (sz * p.dy1 - sy * p.dz1) / (p.r1 * p.r1);
                let cr2 = K * (sy * p.dz2 - sz * p.dy2) / (p.r2 * p.r2);
                cy1 += 2.0 * p.dy1 * cr1 * p.q1 / p.r1;
                cz1 += 2.0 * p.dz1 * cr1 * p.q1 / p.r1;
                cy2 += 2.0 * p.dy2 * cr2 * p.q2 / p.r2;
                cz2 += 2.0 * p.dz2 * cr2 * p.q2 / p.r2;
                let center = 0.5 * Vec3::new(0.0, cy1 + cy2, cz1 + cz2);
                cot.bound_left[target.panel] += center;
                cot.bound_right[target.panel] += center;
                let mut left = Vec3::new(0.0, -cy1, -cz1);
                let mut right = Vec3::new(0.0, -cy2, -cz2);
                let cc = 2.0 * p.core.powi(3) * (cr1 / p.r1 + cr2 / p.r2);
                if p.chord_core {
                    cot.strip_chord[j] += 0.01 * cc;
                } else if p.width > 0.0 {
                    let w = 0.1 * cc / p.width
                        * Vec3::new(0.0, p.right.y - p.left.y, p.right.z - p.left.z);
                    right += w;
                    left -= w;
                }
                if image == 1 {
                    cot.bound_right[source.panel] += reflect(left);
                    cot.bound_left[source.panel] += reflect(right);
                } else {
                    cot.bound_left[source.panel] += left;
                    cot.bound_right[source.panel] += right;
                }
            }
        }
    }
    let mut result = VecX::zeros(s.mesh().vortex_count());
    for (i, strip) in s.mesh().strips().iter().enumerate() {
        for j in strip.first_vortex_index..strip.first_vortex_index + strip.num_vortices {
            result[j] = gamma[i];
        }
    }
    (result, cot)
}
