/*
    Project Pegasus
    Version 0.2.0

    Copyright 2024-2026, Christopher A. Lupp

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

#ifndef PEGASUS_VISCOUS_BORA_POLAR_H_
#define PEGASUS_VISCOUS_BORA_POLAR_H_

#include "pegasus/core/types.h"
#include "pegasus/viscous/airfoil_solver.h"
#include "pegasus/viscous/coupling_station.h"

#include <memory>
#include <string>
#include <vector>

#include <Eigen/Core>

namespace pegasus::viscous {

/// @brief Network size of the Bora (NeuralFoil) surrogate. Mirrors
/// bora::ModelSize; larger networks are more accurate and slower.
enum class BoraModelSize {
    kXXSmall,
    kXSmall,
    kSmall,
    kMedium,
    kLarge,
    kXLarge,  ///< Default in upstream NeuralFoil.
    kXXLarge,
    kXXXLarge,
};

/// @brief Number of Kulfan (CST) shape parameters in Bora's shape space:
/// 8 upper weights, 8 lower weights, leading-edge weight, TE thickness.
inline constexpr int kBoraNumShapeParams = 18;

/// @brief Settings for the Bora airfoil-solver backend, shared by every
/// coupling station.
struct BoraOptions {
    /// Neural network size (speed/accuracy tradeoff).
    BoraModelSize model_size = BoraModelSize::kXLarge;

    /// Critical amplification factor for the transition criterion.
    Scalar n_crit = 9.0;

    /// Forced transition location on the upper surface (x/c); 1.0 = free.
    Scalar xtr_upper = 1.0;

    /// Forced transition location on the lower surface (x/c); 1.0 = free.
    Scalar xtr_lower = 1.0;

    /// Directory containing the .h5 network weight files. Empty uses Bora's
    /// default resolution (BORA_DATA_DIR environment variable, then the
    /// install/source data directories baked in at Bora build time).
    std::string data_dir{};

    /// Minimum analysis confidence in [0, 1] for a solve to count as
    /// converged. The surrogate always produces an output; confidence
    /// measures proximity to its training distribution. 0 accepts every
    /// evaluation.
    Scalar min_confidence = 0.0;
};

/// @brief AirfoilSolver backend wrapping Bora, a C++ port of the NeuralFoil
/// neural surrogate of XFoil.
///
/// Construction fits Kulfan (CST) parameters to the station's airfoil
/// coordinates (any closed-loop orientation; converted to Selig order and
/// chord-normalized internally). Each solve is a single network evaluation
/// with analytic gradients — there is no boundary-layer iteration and no
/// warm-start state, so solves are fast (~ms) and never diverge.
///
/// Limitations relative to the Scirocco backend:
///   - Mach is ignored: the surrogate is trained on incompressible XFoil
///     data (valid for low subsonic sections).
///   - Shape derivatives (GetDclDx / GetDcdDx) are 18×1 in Kulfan space
///     (see kBoraNumShapeParams) rather than per-panel-node coordinates.
class BoraPolarSolver final : public AirfoilSolver {
public:
    /// @brief Construct from a 2×N airfoil node coordinate matrix, loading
    /// a private network model. To share one model across many stations,
    /// use CreateBoraStationSolvers instead.
    /// @throws std::runtime_error if the weight files cannot be loaded or
    /// the Kulfan fit is singular.
    BoraPolarSolver(const Eigen::MatrixXd& coords, const BoraOptions& options);

    ~BoraPolarSolver() override;
    BoraPolarSolver(BoraPolarSolver&&) noexcept;
    BoraPolarSolver& operator=(BoraPolarSolver&&) noexcept;
    BoraPolarSolver(const BoraPolarSolver&) = delete;
    BoraPolarSolver& operator=(const BoraPolarSolver&) = delete;

    /// @brief Evaluate the surrogate at fixed α (rad) and Reynolds. Mach is
    /// ignored (see class docs). Returns false only when the analysis
    /// confidence falls below BoraOptions::min_confidence.
    bool SolveFixedAlpha(Scalar alpha_rad, Scalar reynolds,
                         Scalar mach) override;

    /// @brief Angle of attack of the most recent solve, radians.
    [[nodiscard]] Scalar GetAlphaRadians() const override;

    /// @brief Convergence (confidence) flag of the most recent solve.
    [[nodiscard]] bool IsConverged() const override;

    /// @brief Section lift coefficient of the most recent solve.
    [[nodiscard]] Scalar GetCl() const override;
    /// @brief Section drag coefficient of the most recent solve.
    [[nodiscard]] Scalar GetCd() const override;
    /// @brief Section quarter-chord moment coefficient of the most recent
    /// solve.
    [[nodiscard]] Scalar GetCm() const override;

    /// @brief dCl/dα in units of per-radian (analytic).
    [[nodiscard]] Scalar GetDclDalphaPerRadian() const override;
    /// @brief dCd/dα in units of per-radian (analytic).
    [[nodiscard]] Scalar GetDcdDalphaPerRadian() const override;

    /// @brief dCl / d(Kulfan params), an 18×1 matrix ordered as
    /// [upper 0..7, lower 0..7, LE weight, TE thickness].
    [[nodiscard]] const Eigen::MatrixXd& GetDclDx() const override;
    /// @brief dCd / d(Kulfan params) — same layout as GetDclDx.
    [[nodiscard]] const Eigen::MatrixXd& GetDcdDx() const override;

    /// @brief Network analysis confidence in [0, 1] of the most recent
    /// solve (proximity of the operating point to Bora's training
    /// distribution).
    [[nodiscard]] Scalar GetAnalysisConfidence() const;

private:
    BoraPolarSolver();

    friend std::vector<std::unique_ptr<AirfoilSolver>> CreateBoraStationSolvers(
        const std::vector<CouplingStation>& stations,
        const BoraOptions& options);

    struct Impl;
    std::unique_ptr<Impl> impl_;
};

/// @brief Build one BoraPolarSolver per coupling station, sharing a single
/// neural network model across all stations (the model is immutable and
/// thread-safe, so per-station solves may run concurrently).
/// @throws std::runtime_error if the weight files cannot be loaded or a
/// station's Kulfan fit is singular.
std::vector<std::unique_ptr<AirfoilSolver>> CreateBoraStationSolvers(
    const std::vector<CouplingStation>& stations, const BoraOptions& options);

}  // namespace pegasus::viscous

#endif  // PEGASUS_VISCOUS_BORA_POLAR_H_
