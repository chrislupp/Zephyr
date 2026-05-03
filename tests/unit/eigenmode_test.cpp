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

#include <cmath>

#include <gtest/gtest.h>

#include "vlm/dynamics/eigenmode_analysis.h"
#include "vlm/dynamics/mass_properties.h"
#include "vlm/geometry/aircraft.h"
#include "vlm/geometry/section.h"
#include "vlm/geometry/surface.h"
#include "vlm/solver/vortex_lattice_solver.h"

namespace vlm {
namespace {

constexpr double kTestTolerance = 1.0e-6;

// ============================================================================
// MassProperties Tests
// ============================================================================

class MassPropertiesTest : public ::testing::Test {};

TEST_F(MassPropertiesTest, DefaultIsInvalid) {
  MassProperties props;
  EXPECT_FALSE(props.IsValid());
}

TEST_F(MassPropertiesTest, ValidWithPositiveValues) {
  MassProperties props;
  props.mass = 1000.0;
  props.Ixx = 1000.0;
  props.Iyy = 2000.0;
  props.Izz = 3000.0;

  EXPECT_TRUE(props.IsValid());
}

TEST_F(MassPropertiesTest, InertiaTensorSymmetric) {
  MassProperties props;
  props.Ixx = 100.0;
  props.Iyy = 200.0;
  props.Izz = 300.0;
  props.Ixy = 10.0;
  props.Iyz = 20.0;
  props.Ixz = 30.0;

  Mat3 I = props.GetInertiaTensor();

  // Check symmetry
  EXPECT_NEAR(I(0, 1), I(1, 0), kTestTolerance);
  EXPECT_NEAR(I(1, 2), I(2, 1), kTestTolerance);
  EXPECT_NEAR(I(0, 2), I(2, 0), kTestTolerance);

  // Check diagonal
  EXPECT_NEAR(I(0, 0), 100.0, kTestTolerance);
  EXPECT_NEAR(I(1, 1), 200.0, kTestTolerance);
  EXPECT_NEAR(I(2, 2), 300.0, kTestTolerance);
}

TEST_F(MassPropertiesTest, InverseInertiaCorrect) {
  MassProperties props;
  props.mass = 1000.0;
  props.Ixx = 100.0;
  props.Iyy = 200.0;
  props.Izz = 300.0;

  Mat3 I = props.GetInertiaTensor();
  Mat3 I_inv = props.GetInverseInertiaTensor();
  Mat3 identity = I * I_inv;

  // Should be identity matrix
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      Scalar expected = (i == j) ? 1.0 : 0.0;
      EXPECT_NEAR(identity(i, j), expected, kTestTolerance);
    }
  }
}

// ============================================================================
// FlightState Tests
// ============================================================================

TEST(FlightStateTest, BodyVelocityAtZeroAngles) {
  FlightState state;
  state.velocity = 100.0;
  state.alpha = 0.0;
  state.beta = 0.0;

  Vec3 v_body = state.GetBodyVelocity();

  EXPECT_NEAR(v_body.x(), 100.0, kTestTolerance);
  EXPECT_NEAR(v_body.y(), 0.0, kTestTolerance);
  EXPECT_NEAR(v_body.z(), 0.0, kTestTolerance);
}

TEST(FlightStateTest, BodyVelocityWithAlpha) {
  FlightState state;
  state.velocity = 100.0;
  state.alpha = 5.0 * kDegreesToRadians;
  state.beta = 0.0;

  Vec3 v_body = state.GetBodyVelocity();

  // At positive alpha: u = V*cos(alpha), w = V*sin(alpha)
  EXPECT_NEAR(v_body.x(), 100.0 * std::cos(state.alpha), kTestTolerance);
  EXPECT_NEAR(v_body.y(), 0.0, kTestTolerance);
  EXPECT_NEAR(v_body.z(), 100.0 * std::sin(state.alpha), kTestTolerance);
}

TEST(FlightStateTest, DynamicPressure) {
  FlightState state;
  state.velocity = 50.0;
  state.density = 1.225;

  Scalar q = state.GetDynamicPressure();
  Scalar expected = 0.5 * 1.225 * 50.0 * 50.0;

  EXPECT_NEAR(q, expected, kTestTolerance);
}

// ============================================================================
// EigenmodeAnalysis Tests
// ============================================================================

class EigenmodeAnalysisTest : public ::testing::Test {
 protected:
  // Create typical light aircraft mass properties
  MassProperties CreateTypicalMassProperties() {
    MassProperties props;
    props.mass = 1000.0;     // kg
    props.Ixx = 1000.0;      // kg*m^2 (roll inertia)
    props.Iyy = 2000.0;      // kg*m^2 (pitch inertia)
    props.Izz = 3000.0;      // kg*m^2 (yaw inertia)
    props.Ixy = 0.0;
    props.Iyz = 0.0;
    props.Ixz = 100.0;       // Typical product of inertia
    return props;
  }

  // Create typical flight state
  FlightState CreateTypicalFlightState() {
    FlightState state;
    state.velocity = 50.0;       // m/s
    state.alpha = 5.0 * kDegreesToRadians;
    state.beta = 0.0;
    state.density = 1.225;       // kg/m^3 (sea level)
    state.gravity = 9.81;        // m/s^2
    state.phi = 0.0;
    state.theta = state.alpha;   // Level flight
    return state;
  }

  // Create typical stability derivatives for a stable aircraft
  StabilityDerivatives CreateTypicalDerivatives() {
    StabilityDerivatives derivs;
    derivs.is_valid = true;

    derivs.reference_area = 16.0;   // m^2
    derivs.reference_chord = 1.5;   // m
    derivs.reference_span = 10.0;   // m

    // Longitudinal derivatives
    derivs.CL_alpha = 5.0;          // Lift curve slope (per rad)
    derivs.CD_alpha = 0.1;          // Drag due to alpha
    derivs.Cm_alpha = -1.0;         // Static stability (negative = stable)
    derivs.Cm_q = -15.0;            // Pitch damping (negative = stable)
    derivs.CL_q = 5.0;              // Lift due to pitch rate

    // Lateral-directional derivatives
    derivs.CY_beta = -0.5;          // Side force due to sideslip
    derivs.Cl_beta = -0.1;          // Dihedral effect (negative = stable)
    derivs.Cn_beta = 0.1;           // Directional stability (positive = stable)
    derivs.Cl_p = -0.5;             // Roll damping (negative = stable)
    derivs.Cn_r = -0.1;             // Yaw damping (negative = stable)
    derivs.Cl_r = 0.1;              // Roll due to yaw rate
    derivs.Cn_p = -0.05;            // Yaw due to roll rate
    derivs.CY_p = 0.0;
    derivs.CY_r = 0.2;

    return derivs;
  }
};

TEST_F(EigenmodeAnalysisTest, LongitudinalAnalysisProducesModes) {
  MassProperties mass = CreateTypicalMassProperties();
  FlightState flight = CreateTypicalFlightState();
  StabilityDerivatives derivs = CreateTypicalDerivatives();

  Aircraft aircraft("Test");
  aircraft.SetReferenceArea(derivs.reference_area);
  aircraft.SetReferenceChord(derivs.reference_chord);
  aircraft.SetReferenceSpan(derivs.reference_span);

  EigenmodeAnalysis analysis;
  EigenmodeResults results =
      analysis.AnalyzeLongitudinal(aircraft, mass, derivs, flight);

  EXPECT_TRUE(results.is_valid);
  EXPECT_EQ(results.num_states, 4);
  EXPECT_GE(results.modes.size(), 1);
}

TEST_F(EigenmodeAnalysisTest, LateralDirectionalAnalysisProducesModes) {
  MassProperties mass = CreateTypicalMassProperties();
  FlightState flight = CreateTypicalFlightState();
  StabilityDerivatives derivs = CreateTypicalDerivatives();

  Aircraft aircraft("Test");
  aircraft.SetReferenceArea(derivs.reference_area);
  aircraft.SetReferenceChord(derivs.reference_chord);
  aircraft.SetReferenceSpan(derivs.reference_span);

  EigenmodeAnalysis analysis;
  EigenmodeResults results =
      analysis.AnalyzeLateralDirectional(aircraft, mass, derivs, flight);

  EXPECT_TRUE(results.is_valid);
  EXPECT_EQ(results.num_states, 4);
  EXPECT_GE(results.modes.size(), 1);
}

TEST_F(EigenmodeAnalysisTest, StableAircraftHasNegativeDamping) {
  MassProperties mass = CreateTypicalMassProperties();
  FlightState flight = CreateTypicalFlightState();
  StabilityDerivatives derivs = CreateTypicalDerivatives();

  Aircraft aircraft("Test");
  aircraft.SetReferenceArea(derivs.reference_area);
  aircraft.SetReferenceChord(derivs.reference_chord);
  aircraft.SetReferenceSpan(derivs.reference_span);

  EigenmodeAnalysis analysis;
  EigenmodeResults results =
      analysis.AnalyzeLongitudinal(aircraft, mass, derivs, flight);

  EXPECT_TRUE(results.is_valid);

  // All modes should have negative real parts (stable)
  for (const auto& mode : results.modes) {
    EXPECT_LE(mode.GetDamping(), 0.0)
        << "Mode with eigenvalue " << mode.eigenvalue << " is unstable";
  }
}

TEST_F(EigenmodeAnalysisTest, ShortPeriodHigherFrequencyThanPhugoid) {
  MassProperties mass = CreateTypicalMassProperties();
  FlightState flight = CreateTypicalFlightState();
  StabilityDerivatives derivs = CreateTypicalDerivatives();

  Aircraft aircraft("Test");
  aircraft.SetReferenceArea(derivs.reference_area);
  aircraft.SetReferenceChord(derivs.reference_chord);
  aircraft.SetReferenceSpan(derivs.reference_span);

  EigenmodeAnalysis analysis;
  EigenmodeResults results =
      analysis.AnalyzeLongitudinal(aircraft, mass, derivs, flight);

  const Eigenmode* sp = results.GetShortPeriod();
  const Eigenmode* ph = results.GetPhugoid();

  if (sp && ph) {
    EXPECT_GT(sp->GetFrequency(), ph->GetFrequency());
  }
}

TEST_F(EigenmodeAnalysisTest, EigenmodeGettersWork) {
  Eigenmode mode;
  mode.eigenvalue = ComplexScalar(-0.5, 2.0);

  EXPECT_NEAR(mode.GetDamping(), -0.5, kTestTolerance);
  EXPECT_NEAR(mode.GetFrequency(), 2.0, kTestTolerance);
  EXPECT_TRUE(mode.IsStable());
  EXPECT_TRUE(mode.IsOscillatory());

  Scalar period = mode.GetPeriod();
  EXPECT_NEAR(period, kPi, 0.01);  // Period = 2*pi/omega = 2*pi/2 = pi
}

TEST_F(EigenmodeAnalysisTest, NonOscillatoryModeHasZeroPeriod) {
  Eigenmode mode;
  mode.eigenvalue = ComplexScalar(-1.0, 0.0);

  EXPECT_FALSE(mode.IsOscillatory());
  EXPECT_NEAR(mode.GetPeriod(), 0.0, kTestTolerance);
}

TEST_F(EigenmodeAnalysisTest, DampingRatioComputation) {
  Eigenmode mode;
  // For critically damped: zeta = 1, sigma = -omega_n, omega = 0
  // For typical: sigma = -zeta*omega_n, omega = omega_n*sqrt(1-zeta^2)

  // Under-damped case: zeta = 0.5, omega_n = 2
  Scalar zeta = 0.5;
  Scalar omega_n = 2.0;
  Scalar sigma = -zeta * omega_n;
  Scalar omega = omega_n * std::sqrt(1 - zeta * zeta);

  mode.eigenvalue = ComplexScalar(sigma, omega);

  EXPECT_NEAR(mode.GetDampingRatio(), zeta, 0.01);
}

TEST_F(EigenmodeAnalysisTest, InvalidMassPropertiesFails) {
  MassProperties mass;  // Invalid - all zeros
  FlightState flight = CreateTypicalFlightState();
  StabilityDerivatives derivs = CreateTypicalDerivatives();

  Aircraft aircraft("Test");

  EigenmodeAnalysis analysis;
  EigenmodeResults results =
      analysis.AnalyzeLongitudinal(aircraft, mass, derivs, flight);

  EXPECT_FALSE(results.is_valid);
}

TEST_F(EigenmodeAnalysisTest, InvalidDerivativesFails) {
  MassProperties mass = CreateTypicalMassProperties();
  FlightState flight = CreateTypicalFlightState();
  StabilityDerivatives derivs;  // Invalid - is_valid = false

  Aircraft aircraft("Test");

  EigenmodeAnalysis analysis;
  EigenmodeResults results =
      analysis.AnalyzeLongitudinal(aircraft, mass, derivs, flight);

  EXPECT_FALSE(results.is_valid);
}

TEST_F(EigenmodeAnalysisTest, SystemMatrixHasCorrectSize) {
  MassProperties mass = CreateTypicalMassProperties();
  FlightState flight = CreateTypicalFlightState();
  StabilityDerivatives derivs = CreateTypicalDerivatives();

  Aircraft aircraft("Test");
  aircraft.SetReferenceArea(derivs.reference_area);
  aircraft.SetReferenceChord(derivs.reference_chord);
  aircraft.SetReferenceSpan(derivs.reference_span);

  EigenmodeAnalysis analysis;
  EigenmodeResults results =
      analysis.AnalyzeLongitudinal(aircraft, mass, derivs, flight);

  EXPECT_EQ(results.system_matrix.rows(), 4);
  EXPECT_EQ(results.system_matrix.cols(), 4);
}

TEST_F(EigenmodeAnalysisTest, FullAnalysisProducesMoreModes) {
  MassProperties mass = CreateTypicalMassProperties();
  FlightState flight = CreateTypicalFlightState();
  StabilityDerivatives derivs = CreateTypicalDerivatives();

  Aircraft aircraft("Test");
  aircraft.SetReferenceArea(derivs.reference_area);
  aircraft.SetReferenceChord(derivs.reference_chord);
  aircraft.SetReferenceSpan(derivs.reference_span);

  EigenmodeAnalysis analysis;
  EigenmodeResults full_results = analysis.Analyze(aircraft, mass, derivs, flight);
  EigenmodeResults lon_results =
      analysis.AnalyzeLongitudinal(aircraft, mass, derivs, flight);

  EXPECT_TRUE(full_results.is_valid);
  EXPECT_GT(full_results.num_states, lon_results.num_states);
}

TEST_F(EigenmodeAnalysisTest, ModeTypeNames) {
  EXPECT_EQ(GetModeTypeName(ModeType::kShortPeriod), "Short Period");
  EXPECT_EQ(GetModeTypeName(ModeType::kPhugoid), "Phugoid");
  EXPECT_EQ(GetModeTypeName(ModeType::kDutchRoll), "Dutch Roll");
  EXPECT_EQ(GetModeTypeName(ModeType::kRollSubsidence), "Roll Subsidence");
  EXPECT_EQ(GetModeTypeName(ModeType::kSpiralMode), "Spiral");
}

}  // namespace
}  // namespace vlm
