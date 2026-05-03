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

#ifndef VLM_DYNAMICS_EIGENMODE_ANALYSIS_H_
#define VLM_DYNAMICS_EIGENMODE_ANALYSIS_H_

#include "vlm/aero/stability_derivatives.h"
#include "vlm/core/types.h"
#include "vlm/dynamics/mass_properties.h"
#include "vlm/geometry/aircraft.h"

#include <complex>
#include <string>
#include <vector>

namespace vlm {

/// Complex eigenvalue type
using ComplexScalar = std::complex<Scalar>;

/// Classification of dynamic mode
enum class ModeType {
    kUnknown,
    kShortPeriod,     ///< High-frequency pitch oscillation
    kPhugoid,         ///< Low-frequency speed/altitude oscillation
    kDutchRoll,       ///< Coupled yaw/roll oscillation
    kRollSubsidence,  ///< Pure roll damping (real eigenvalue)
    kSpiralMode,      ///< Slow heading divergence/convergence
    kOther            ///< Unclassified mode
};

/// Get string name for mode type
inline std::string GetModeTypeName(ModeType type) {
    switch (type) {
        case ModeType::kShortPeriod:
            return "Short Period";
        case ModeType::kPhugoid:
            return "Phugoid";
        case ModeType::kDutchRoll:
            return "Dutch Roll";
        case ModeType::kRollSubsidence:
            return "Roll Subsidence";
        case ModeType::kSpiralMode:
            return "Spiral";
        default:
            return "Unknown";
    }
}

/// Single eigenmode result
struct Eigenmode {
    ComplexScalar eigenvalue{0.0,
                             0.0};  ///< Complex eigenvalue (sigma + i*omega)

    /// Eigenvector components (complex, normalized)
    /// Order: [u, v, w, p, q, r, phi, theta, psi, x, y, z]
    Eigen::Matrix<ComplexScalar, Eigen::Dynamic, 1> eigenvector;

    ModeType type{ModeType::kUnknown};  ///< Classified mode type

    /// Real part (damping): negative = stable
    Scalar GetDamping() const { return eigenvalue.real(); }

    /// Imaginary part (frequency in rad/s)
    Scalar GetFrequency() const { return std::abs(eigenvalue.imag()); }

    /// Period in seconds (0 if non-oscillatory)
    Scalar GetPeriod() const {
        Scalar omega = GetFrequency();
        return omega > 0 ? 2.0 * kPi / omega : 0.0;
    }

    /// Time to half/double amplitude (positive = time to half)
    Scalar GetTimeConstant() const {
        Scalar sigma = GetDamping();
        return sigma != 0 ? -std::log(2.0) / sigma : 0.0;
    }

    /// Damping ratio (for oscillatory modes)
    Scalar GetDampingRatio() const {
        Scalar sigma = GetDamping();
        Scalar omega = GetFrequency();
        Scalar omega_n = std::sqrt(sigma * sigma + omega * omega);
        return omega_n > 0 ? -sigma / omega_n : 0.0;
    }

    /// Check if mode is stable (real part < 0)
    bool IsStable() const { return eigenvalue.real() < 0; }

    /// Check if mode is oscillatory (has imaginary part)
    bool IsOscillatory() const {
        return std::abs(eigenvalue.imag()) > kEpsilon;
    }
};

/// Complete eigenmode analysis results
struct EigenmodeResults {
    bool is_valid{false};

    /// All computed eigenmodes
    std::vector<Eigenmode> modes;

    /// System matrix A (dx/dt = A*x)
    MatX system_matrix;

    /// Number of states in the system
    Index num_states{0};

    /// Reference quantities used
    Scalar reference_velocity{0.0};
    Scalar reference_span{0.0};
    Scalar reference_chord{0.0};

    /// Find modes by type
    std::vector<const Eigenmode*> FindModesByType(ModeType type) const {
        std::vector<const Eigenmode*> result;
        for (const auto& mode : modes) {
            if (mode.type == type) {
                result.push_back(&mode);
            }
        }
        return result;
    }

    /// Get short period mode (if exists)
    const Eigenmode* GetShortPeriod() const {
        auto found = FindModesByType(ModeType::kShortPeriod);
        return found.empty() ? nullptr : found[0];
    }

    /// Get phugoid mode (if exists)
    const Eigenmode* GetPhugoid() const {
        auto found = FindModesByType(ModeType::kPhugoid);
        return found.empty() ? nullptr : found[0];
    }

    /// Get Dutch roll mode (if exists)
    const Eigenmode* GetDutchRoll() const {
        auto found = FindModesByType(ModeType::kDutchRoll);
        return found.empty() ? nullptr : found[0];
    }

    /// Check if all modes are stable
    bool IsStable() const {
        for (const auto& mode : modes) {
            if (!mode.IsStable())
                return false;
        }
        return true;
    }
};

/// Eigenmode analysis for aircraft dynamic stability
///
/// Computes the linearized equations of motion and solves for eigenvalues
/// to determine dynamic stability characteristics. The analysis uses:
///
/// - 6 DOF rigid body equations
/// - Aerodynamic derivatives from VLM
/// - Mass and inertia properties
///
/// State vector (12 states):
/// [u, v, w, p, q, r, phi, theta, psi, x, y, z]
/// where u,v,w are body-axis velocity perturbations,
/// p,q,r are body-axis angular rates,
/// phi,theta,psi are Euler angles,
/// x,y,z are position.
///
/// For lateral-directional analysis, a reduced 4-state model can be used:
/// [v, p, r, phi] for Dutch roll, spiral, and roll modes.
///
/// For longitudinal analysis, a reduced 4-state model can be used:
/// [u, w, q, theta] for short period and phugoid modes.
class EigenmodeAnalysis {
public:
    EigenmodeAnalysis() = default;

    /// Perform full 12-state eigenmode analysis
    ///
    /// @param aircraft Aircraft geometry
    /// @param mass_props Mass and inertia properties
    /// @param derivs Stability derivatives at trim condition
    /// @param flight Flight state at trim
    /// @return Complete eigenmode results
    EigenmodeResults Analyze(const Aircraft& aircraft,
                             const MassProperties& mass_props,
                             const StabilityDerivatives& derivs,
                             const FlightState& flight);

    /// Perform longitudinal-only analysis (4-state)
    ///
    /// Analyzes only the symmetric (longitudinal) modes:
    /// - Short period
    /// - Phugoid
    ///
    /// @param aircraft Aircraft geometry
    /// @param mass_props Mass and inertia properties
    /// @param derivs Stability derivatives at trim condition
    /// @param flight Flight state at trim
    /// @return Longitudinal eigenmode results
    EigenmodeResults AnalyzeLongitudinal(const Aircraft& aircraft,
                                         const MassProperties& mass_props,
                                         const StabilityDerivatives& derivs,
                                         const FlightState& flight);

    /// Perform lateral-directional-only analysis (4-state)
    ///
    /// Analyzes only the asymmetric (lateral-directional) modes:
    /// - Dutch roll
    /// - Roll subsidence
    /// - Spiral mode
    ///
    /// @param aircraft Aircraft geometry
    /// @param mass_props Mass and inertia properties
    /// @param derivs Stability derivatives at trim condition
    /// @param flight Flight state at trim
    /// @return Lateral-directional eigenmode results
    EigenmodeResults AnalyzeLateralDirectional(
        const Aircraft& aircraft, const MassProperties& mass_props,
        const StabilityDerivatives& derivs, const FlightState& flight);

    /// Set apparent mass effects
    void SetApparentMass(const ApparentMass& apparent_mass) {
        apparent_mass_ = apparent_mass;
    }

private:
    /// Build the full 12x12 system matrix
    MatX BuildFullSystemMatrix(const Aircraft& aircraft,
                               const MassProperties& mass_props,
                               const StabilityDerivatives& derivs,
                               const FlightState& flight) const;

    /// Build longitudinal 4x4 system matrix [u, w, q, theta]
    Mat4 BuildLongitudinalMatrix(const Aircraft& aircraft,
                                 const MassProperties& mass_props,
                                 const StabilityDerivatives& derivs,
                                 const FlightState& flight) const;

    /// Build lateral-directional 4x4 system matrix [v, p, r, phi]
    Mat4 BuildLateralDirectionalMatrix(const Aircraft& aircraft,
                                       const MassProperties& mass_props,
                                       const StabilityDerivatives& derivs,
                                       const FlightState& flight) const;

    /// Compute eigenvalues and eigenvectors
    void ComputeEigenvalues(const MatX& system_matrix,
                            std::vector<Eigenmode>& modes) const;

    /// Classify mode type based on eigenvector content and frequency
    ModeType ClassifyMode(const Eigenmode& mode, bool is_longitudinal) const;

    ApparentMass apparent_mass_;
};

}  // namespace vlm

#endif  // VLM_DYNAMICS_EIGENMODE_ANALYSIS_H_
