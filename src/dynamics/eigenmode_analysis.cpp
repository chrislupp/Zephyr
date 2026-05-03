/*
    Project VLM
    Version 0.0.0

    Copyright 2024, Christopher A. Lupp

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published
    by the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include "vlm/dynamics/eigenmode_analysis.h"

#include <algorithm>
#include <cmath>

#include <Eigen/Eigenvalues>

namespace vlm {

namespace {

/// Compute rotation matrix from Euler angles (3-2-1 sequence: yaw-pitch-roll)
Mat3 EulerToRotationMatrix(Scalar phi, Scalar theta,
                           [[maybe_unused]] Scalar psi) {
    Scalar cp = std::cos(phi);
    Scalar sp = std::sin(phi);
    Scalar ct = std::cos(theta);
    Scalar st = std::sin(theta);

    // Rotation matrix (body to inertial)
    Mat3 R;
    R << ct, sp * st, cp * st, 0, cp, -sp, -st, sp * ct, cp * ct;
    return R;
}

/// Compute rate transformation matrix (body rates to Euler rates)
Mat3 BodyRatesToEulerRates(Scalar phi, Scalar theta) {
    Scalar cp = std::cos(phi);
    Scalar sp = std::sin(phi);
    Scalar ct = std::cos(theta);
    Scalar tt = std::tan(theta);

    Mat3 L;
    L << 1, sp * tt, cp * tt, 0, cp, -sp, 0, sp / ct, cp / ct;
    return L;
}

}  // namespace

EigenmodeResults EigenmodeAnalysis::Analyze(const Aircraft& aircraft,
                                            const MassProperties& mass_props,
                                            const StabilityDerivatives& derivs,
                                            const FlightState& flight) {
    EigenmodeResults results;

    if (!mass_props.IsValid() || !derivs.is_valid) {
        return results;
    }

    // Build system matrix
    results.system_matrix =
        BuildFullSystemMatrix(aircraft, mass_props, derivs, flight);
    results.num_states = results.system_matrix.rows();
    results.reference_velocity = flight.velocity;
    results.reference_span = derivs.reference_span;
    results.reference_chord = derivs.reference_chord;

    // Compute eigenvalues
    ComputeEigenvalues(results.system_matrix, results.modes);

    // Classify modes
    for (auto& mode : results.modes) {
        mode.type = ClassifyMode(mode, false);
    }

    results.is_valid = true;
    return results;
}

EigenmodeResults EigenmodeAnalysis::AnalyzeLongitudinal(
    const Aircraft& aircraft, const MassProperties& mass_props,
    const StabilityDerivatives& derivs, const FlightState& flight) {
    EigenmodeResults results;

    if (!mass_props.IsValid() || !derivs.is_valid) {
        return results;
    }

    // Build 4x4 longitudinal system matrix
    Mat4 A_lon = BuildLongitudinalMatrix(aircraft, mass_props, derivs, flight);

    results.system_matrix = A_lon;
    results.num_states = 4;
    results.reference_velocity = flight.velocity;
    results.reference_span = derivs.reference_span;
    results.reference_chord = derivs.reference_chord;

    // Compute eigenvalues
    ComputeEigenvalues(results.system_matrix, results.modes);

    // Classify modes
    for (auto& mode : results.modes) {
        mode.type = ClassifyMode(mode, true);
    }

    results.is_valid = true;
    return results;
}

EigenmodeResults EigenmodeAnalysis::AnalyzeLateralDirectional(
    const Aircraft& aircraft, const MassProperties& mass_props,
    const StabilityDerivatives& derivs, const FlightState& flight) {
    EigenmodeResults results;

    if (!mass_props.IsValid() || !derivs.is_valid) {
        return results;
    }

    // Build 4x4 lateral-directional system matrix
    Mat4 A_lat =
        BuildLateralDirectionalMatrix(aircraft, mass_props, derivs, flight);

    results.system_matrix = A_lat;
    results.num_states = 4;
    results.reference_velocity = flight.velocity;
    results.reference_span = derivs.reference_span;
    results.reference_chord = derivs.reference_chord;

    // Compute eigenvalues
    ComputeEigenvalues(results.system_matrix, results.modes);

    // Classify modes (lateral modes)
    for (auto& mode : results.modes) {
        mode.type = ClassifyMode(mode, false);
    }

    results.is_valid = true;
    return results;
}

MatX EigenmodeAnalysis::BuildFullSystemMatrix(
    const Aircraft& aircraft, const MassProperties& mass_props,
    const StabilityDerivatives& derivs, const FlightState& flight) const {
    // Full 12-state system: [u, v, w, p, q, r, phi, theta, psi, x, y, z]
    // Simplified to 9-state for practical use: [u, v, w, p, q, r, phi, theta,
    // psi]
    constexpr Index kNumStates = 9;
    MatX A = MatX::Zero(kNumStates, kNumStates);

    Scalar V = flight.velocity;
    Scalar g = flight.gravity;
    Scalar m = mass_props.mass;
    Scalar rho = flight.density;

    Scalar S_ref = aircraft.GetReferenceArea();
    Scalar c_ref = aircraft.GetReferenceChord();
    Scalar b_ref = aircraft.GetReferenceSpan();

    // Dynamic pressure and dimensional conversion factors
    Scalar qS = 0.5 * rho * V * V * S_ref;
    Scalar qS_mV = qS / (m * V);
    Scalar qSb_Ix = qS * b_ref / mass_props.Ixx;
    Scalar qSc_Iy = qS * c_ref / mass_props.Iyy;
    Scalar qSb_Iz = qS * b_ref / mass_props.Izz;

    // State indices
    constexpr Index kU = 0, kV = 1, kW = 2;
    constexpr Index kP = 3, kQ = 4, kR = 5;
    constexpr Index kPhi = 6, kTheta = 7, kPsi = 8;

    // Longitudinal derivatives (u, w, q, theta)
    // du/dt equation (X-force)
    A(kU, kU) = -qS_mV * derivs.CD_alpha * std::sin(flight.alpha);  // X_u
    A(kU, kW) = qS_mV * (derivs.CL_alpha - derivs.CD_alpha);  // X_w approx
    A(kU, kQ) = 0.0;
    A(kU, kTheta) = -g * std::cos(flight.theta);

    // dw/dt equation (Z-force)
    A(kW, kU) = -qS_mV * derivs.CL_alpha * std::cos(flight.alpha);
    A(kW, kW) = -qS_mV * derivs.CL_alpha;
    A(kW, kQ) = V;  // w-dot includes V*q for pitch rate effect
    A(kW, kTheta) = g * std::sin(flight.theta);

    // dq/dt equation (pitch moment)
    A(kQ, kU) = qSc_Iy * derivs.Cm_alpha * (-std::sin(flight.alpha) / V);
    A(kQ, kW) = qSc_Iy * derivs.Cm_alpha / V;
    A(kQ, kQ) = qSc_Iy * derivs.Cm_q * c_ref / (2.0 * V);

    // dtheta/dt = q (kinematic)
    A(kTheta, kQ) = 1.0;

    // Lateral-directional derivatives (v, p, r, phi)
    // dv/dt equation (Y-force)
    A(kV, kV) = qS_mV * derivs.CY_beta;
    A(kV, kP) = qS_mV * derivs.CY_p * b_ref / (2.0 * V);
    A(kV, kR) = qS_mV * derivs.CY_r * b_ref / (2.0 * V) - V;  // includes -V*r
    A(kV, kPhi) = g * std::cos(flight.theta);

    // dp/dt equation (roll moment)
    A(kP, kV) = qSb_Ix * derivs.Cl_beta / V;
    A(kP, kP) = qSb_Ix * derivs.Cl_p * b_ref / (2.0 * V);
    A(kP, kR) = qSb_Ix * derivs.Cl_r * b_ref / (2.0 * V);

    // dr/dt equation (yaw moment)
    A(kR, kV) = qSb_Iz * derivs.Cn_beta / V;
    A(kR, kP) = qSb_Iz * derivs.Cn_p * b_ref / (2.0 * V);
    A(kR, kR) = qSb_Iz * derivs.Cn_r * b_ref / (2.0 * V);

    // dphi/dt = p (kinematic, simplified for small angles)
    A(kPhi, kP) = 1.0;

    // dpsi/dt = r/cos(theta) (kinematic, simplified)
    A(kPsi, kR) = 1.0 / std::cos(flight.theta);

    return A;
}

Mat4 EigenmodeAnalysis::BuildLongitudinalMatrix(
    const Aircraft& aircraft, const MassProperties& mass_props,
    const StabilityDerivatives& derivs, const FlightState& flight) const {
    // 4-state longitudinal: [u, w, q, theta]
    Mat4 A = Mat4::Zero();

    Scalar V = flight.velocity;
    Scalar g = flight.gravity;
    Scalar m = mass_props.mass;
    Scalar rho = flight.density;

    Scalar S_ref = aircraft.GetReferenceArea();
    Scalar c_ref = aircraft.GetReferenceChord();

    Scalar qS = 0.5 * rho * V * V * S_ref;
    Scalar qS_mV = qS / (m * V);
    Scalar qSc_Iy = qS * c_ref / mass_props.Iyy;

    // State indices: u=0, w=1, q=2, theta=3

    // du/dt equation
    A(0, 0) = -qS_mV * 2.0 * derivs.CD_alpha;  // X_u (approx from CD)
    A(0, 1) = qS_mV * derivs.CL_alpha;         // X_w (approx)
    A(0, 2) = 0.0;
    A(0, 3) = -g;

    // dw/dt equation
    A(1, 0) = -qS_mV * derivs.CL_alpha;
    A(1, 1) = -qS_mV * derivs.CL_alpha;
    A(1, 2) = V;
    A(1, 3) = 0.0;

    // dq/dt equation
    A(2, 0) = 0.0;
    A(2, 1) = qSc_Iy * derivs.Cm_alpha / V;
    A(2, 2) = qSc_Iy * derivs.Cm_q * c_ref / (2.0 * V);
    A(2, 3) = 0.0;

    // dtheta/dt = q
    A(3, 0) = 0.0;
    A(3, 1) = 0.0;
    A(3, 2) = 1.0;
    A(3, 3) = 0.0;

    return A;
}

Mat4 EigenmodeAnalysis::BuildLateralDirectionalMatrix(
    const Aircraft& aircraft, const MassProperties& mass_props,
    const StabilityDerivatives& derivs, const FlightState& flight) const {
    // 4-state lateral-directional: [v, p, r, phi]
    Mat4 A = Mat4::Zero();

    Scalar V = flight.velocity;
    Scalar g = flight.gravity;
    Scalar m = mass_props.mass;
    Scalar rho = flight.density;

    Scalar S_ref = aircraft.GetReferenceArea();
    Scalar b_ref = aircraft.GetReferenceSpan();

    Scalar qS = 0.5 * rho * V * V * S_ref;
    Scalar qS_mV = qS / (m * V);
    Scalar qSb_Ix = qS * b_ref / mass_props.Ixx;
    Scalar qSb_Iz = qS * b_ref / mass_props.Izz;

    // State indices: v=0, p=1, r=2, phi=3

    // dv/dt equation (Y-force)
    A(0, 0) = qS_mV * derivs.CY_beta;
    A(0, 1) = qS_mV * derivs.CY_p * b_ref / (2.0 * V);
    A(0, 2) = qS_mV * derivs.CY_r * b_ref / (2.0 * V) - V;
    A(0, 3) = g * std::cos(flight.theta);

    // dp/dt equation (roll moment)
    A(1, 0) = qSb_Ix * derivs.Cl_beta / V;
    A(1, 1) = qSb_Ix * derivs.Cl_p * b_ref / (2.0 * V);
    A(1, 2) = qSb_Ix * derivs.Cl_r * b_ref / (2.0 * V);
    A(1, 3) = 0.0;

    // dr/dt equation (yaw moment)
    A(2, 0) = qSb_Iz * derivs.Cn_beta / V;
    A(2, 1) = qSb_Iz * derivs.Cn_p * b_ref / (2.0 * V);
    A(2, 2) = qSb_Iz * derivs.Cn_r * b_ref / (2.0 * V);
    A(2, 3) = 0.0;

    // dphi/dt = p
    A(3, 0) = 0.0;
    A(3, 1) = 1.0;
    A(3, 2) = 0.0;
    A(3, 3) = 0.0;

    return A;
}

void EigenmodeAnalysis::ComputeEigenvalues(
    const MatX& system_matrix, std::vector<Eigenmode>& modes) const {
    modes.clear();

    Eigen::EigenSolver<MatX> solver(system_matrix);

    auto eigenvalues = solver.eigenvalues();
    auto eigenvectors = solver.eigenvectors();

    Index n = system_matrix.rows();

    // Process eigenvalues, grouping complex conjugate pairs
    std::vector<bool> processed(n, false);

    for (Index i = 0; i < n; ++i) {
        if (processed[i])
            continue;

        Eigenmode mode;
        mode.eigenvalue = eigenvalues(i);
        mode.eigenvector = eigenvectors.col(i);

        // Check for complex conjugate pair
        if (std::abs(eigenvalues(i).imag()) > kEpsilon) {
            // Find conjugate
            for (Index j = i + 1; j < n; ++j) {
                if (!processed[j]) {
                    ComplexScalar diff =
                        eigenvalues(i) - std::conj(eigenvalues(j));
                    if (std::abs(diff) < kEpsilon) {
                        processed[j] = true;
                        break;
                    }
                }
            }
        }

        modes.push_back(mode);
        processed[i] = true;
    }

    // Sort modes by frequency (oscillatory first, then by damping)
    std::sort(modes.begin(), modes.end(),
              [](const Eigenmode& a, const Eigenmode& b) {
                  if (a.IsOscillatory() != b.IsOscillatory()) {
                      return a.IsOscillatory();  // Oscillatory modes first
                  }
                  return a.GetFrequency() >
                         b.GetFrequency();  // Higher frequency first
              });
}

ModeType EigenmodeAnalysis::ClassifyMode(const Eigenmode& mode,
                                         bool is_longitudinal) const {
    Scalar omega = mode.GetFrequency();
    Scalar sigma = mode.GetDamping();
    Scalar zeta = mode.GetDampingRatio();

    if (is_longitudinal) {
        // Longitudinal modes
        if (mode.IsOscillatory()) {
            // Short period: high frequency, well-damped
            // Phugoid: low frequency, lightly damped
            if (omega > 1.0 && zeta > 0.1) {
                return ModeType::kShortPeriod;
            } else if (omega < 1.0) {
                return ModeType::kPhugoid;
            }
        }
    } else {
        // Lateral-directional modes
        if (mode.IsOscillatory()) {
            // Dutch roll: coupled yaw/roll oscillation
            return ModeType::kDutchRoll;
        } else {
            // Non-oscillatory modes
            if (sigma < -0.5) {
                // Fast decay = roll subsidence
                return ModeType::kRollSubsidence;
            } else {
                // Slow mode = spiral
                return ModeType::kSpiralMode;
            }
        }
    }

    return ModeType::kOther;
}

}  // namespace vlm
