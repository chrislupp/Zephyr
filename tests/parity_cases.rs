//! Inventory of the frozen v0.2.0 C++ CTest suite.
//!
//! `tests/fixtures/parity/ctest-manifest.json` is a frozen record of the 488
//! non-Philote cases that made up the C++ suite at tag `v0.2.0` (commit
//! `556b0ea`). It is an *inventory*, not an oracle: it lists which cases
//! existed and what tolerance each ran at, but carries no expected values.
//!
//! This file previously generated one `#[test]` per manifest entry. Those 488
//! tests all dispatched to one of eight category-level smoke checks, so every
//! case in a category ran byte-identical assertions regardless of the C++
//! behaviour it was named after — for example, all 157 `geometry` cases
//! asserted only that a NACA 2412 is valid and has positive camber at x/c=0.4,
//! including the cases named after cubic-spline interpolation. That reported
//! 488 passing parity tests while making eight assertions, and it inflated the
//! line-coverage gate. The generated tests have been removed rather than left
//! to imply coverage that does not exist. The manifest's `rust_test` field
//! records the name of the removed stub for each case.
//!
//! Real behavioural verification lives in the dedicated suites, which assert
//! actual numbers rather than category smoke:
//!
//! - `tests/oracle_parity.rs` — frozen v0.2.0 numerical fixtures
//! - `tests/solver.rs`, `tests/aero.rs`, `tests/gradients.rs`,
//!   `tests/trim.rs`, `tests/dynamics.rs`, `tests/geometry.rs`,
//!   `tests/io.rs`, `tests/visualization.rs`
//!
//! Porting the individual C++ assertions is tracked by [`PORTED`] below. The
//! count is asserted so that it can only move deliberately.

use std::collections::HashSet;

use serde::Deserialize;

#[derive(Debug, Deserialize)]
struct Manifest {
    provenance: Provenance,
    cases: Vec<Case>,
}

#[derive(Debug, Deserialize)]
struct Provenance {
    tag: String,
    commit: String,
    ctest_cases: usize,
    deferred_philote_cases: Vec<String>,
    rust_parity_cases: usize,
}

#[derive(Debug, Deserialize)]
struct Case {
    cpp_test: String,
    category: String,
    tolerance: f64,
    enabled: bool,
}

/// C++ cases whose assertions have been ported to a real Rust test.
///
/// Add an entry only when a Rust test reproduces that case's actual
/// assertions. Membership is verified against the manifest, so a typo or a
/// case that never existed upstream is a hard error.
const PORTED: &[&str] = &[
    // aerodynamic_results_test.cpp -> tests/aero.rs
    "AerodynamicResultsTest.DefaultState",
    "AerodynamicResultsTest.GetAspectRatio",
    "AerodynamicResultsTest.GetAspectRatioSquareWing",
    "AerodynamicResultsTest.GetIdealInducedDrag",
    "AerodynamicResultsTest.GetIdealInducedDragZeroCL",
    "AerodynamicResultsTest.Clear",
    "StripLoadsTest.DefaultState",
    "SurfaceLoadsTest.DefaultState",
    // math_utilities_test.cpp -> tests/math_utilities.rs
    "MathUtilitiesTest.PrandtlGlauertZeroMach",
    "MathUtilitiesTest.PrandtlGlauertSubsonic",
    "MathUtilitiesTest.VortexSegmentAboveCenter",
    "MathUtilitiesTest.VortexSegmentOnAxis",
    "MathUtilitiesTest.VortexSegmentAtEndpoint",
    "MathUtilitiesTest.TrailingVortexAbove",
    "MathUtilitiesTest.TrailingVortexAtStart",
    "MathUtilitiesTest.TrailingVortexDownstream",
    "MathUtilitiesTest.CosineSpacingEndpoints",
    "MathUtilitiesTest.CosineSpacingConcentration",
    "MathUtilitiesTest.UniformSpacing",
    "MathUtilitiesTest.SineSpacing",
    "MathUtilitiesTest.BlendedSpacingUniform",
    "MathUtilitiesTest.BlendedSpacingCosine",
    "MathUtilitiesTest.SpacingDispatcherUniform",
    "MathUtilitiesTest.SpacingDispatcherCosine",
    "MathUtilitiesTest.SpacingDispatcherSine",
    "MathUtilitiesTest.SpacingDispatcherBlend",
    "MathUtilitiesTest.SpacingDispatcherNone",
    "MathUtilitiesTest.RotateAboutXIdentity",
    "MathUtilitiesTest.RotateAboutX90",
    "MathUtilitiesTest.RotateAboutYIdentity",
    "MathUtilitiesTest.RotateAboutY90",
    "MathUtilitiesTest.RotateAboutZIdentity",
    "MathUtilitiesTest.RotateAboutZ90",
    "MathUtilitiesTest.WindToBodyIdentity",
    "MathUtilitiesTest.WindToBodyOrthogonality",
    "MathUtilitiesTest.WindToBodyAlphaOnly",
    "MathUtilitiesTest.BodyToStabilityIdentity",
    "MathUtilitiesTest.BodyToStabilityOrthogonality",
    "MathUtilitiesTest.BodyToStabilityElements",
    // cubic_spline_test.cpp -> tests/cubic_spline.rs
    "CubicSplineTest.LinearDataInterpolation",
    "CubicSplineTest.QuadraticDataWithMatchingBoundaryConditions",
    "CubicSplineTest.SinusoidalData",
    "CubicSplineTest.EvaluateAllConsistency",
    "CubicSplineTest.EndpointInterpolation",
    "CubicSplineTest.SpecifiedSlopeBoundaryCondition",
    "CubicSplineTest.TwoPointSpline",
    "CubicSplineTest.NonUniformSpacing",
    "CubicSplineTest.ExtrapolationBehavior",
    "CubicSplineTest.InvalidInputThrows",
    "CubicSplineTest.EvaluateBeforeBuildThrows",
    "CubicSplineTest.ValidityAndSizeQueries",
    // body_test.cpp -> tests/body.rs
    "BodyTest.DefaultConstruction",
    "BodyTest.NamedConstruction",
    "BodyTest.SetName",
    "BodyTest.AddAndGetSections",
    "BodyTest.NodeCount",
    "BodyTest.NodeSpacing",
    "BodyTest.Translation",
    "BodyTest.Scale",
    "BodyTest.DefaultScale",
    "BodyTest.DefaultTranslation",
    "BodyTest.YDuplicate",
    // mesh_input_test.cpp -> tests/mesh_input.rs
    "MeshInputTest.BasicSolveProducesValidResults",
    "MeshInputTest.ConsistencyWithAircraftPath",
    "MeshInputTest.MeshDimensions",
    // trim_solver_test.cpp -> tests/trim_solver.rs
    "TrimSolverTest.TrimToTargetCL",
    "TrimSolverTest.TrimToZeroCL",
    "TrimSolverTest.TrimWithFixedAlpha",
    "TrimSolverTest.TrimWithFixedBeta",
    "TrimSolverTest.ConvergenceInFewIterations",
    "TrimSolverTest.RatesRemainZeroForLevelFlight",
    "TrimSolverTest.InitialGuessIsUsed",
    "TrimSolverTest.ResidualDecreasesMonotonically",
    "TrimSolverTest.StabilityDerivativesComputed",
    "TrimSolverTest.NoConstraintsFails",
    "TrimSolverTest.UnreachableCLFails",
    "TrimSolverTest.PresetLevelFlight",
    "TrimSolverTest.PresetSteadySideslip",
    "TrimSolverTest.TrimWithPitchMomentConstraint",
    "TrimSolverTest.TrimVerboseOutput",
    "TrimSolverTest.TrimMaxIterationsReached",
    "TrimSolverTest.ConstraintCount",
    // influence_matrix_test.cpp -> tests/influence_matrix.rs
    "BiotSavartTest.VelocityAboveBoundVortex",
    "BiotSavartTest.VelocitySymmetry",
    "BiotSavartTest.CoreRegularization",
    "BiotSavartTest.PrandtlGlauertCorrection",
    "BiotSavartTest.ExcludeBoundVortex",
    "BiotSavartTest.TrailingVortexInDownstreamField",
    "InfluenceMatrixTest.VelocityMatrixDimensions",
    "InfluenceMatrixTest.AICMatrixDimensions",
    "InfluenceMatrixTest.AICDiagonalDominance",
    "InfluenceMatrixTest.AICMatrixNonsingular",
    "InfluenceMatrixTest.MachNumberEffect",
    "InfluenceMatrixTest.InducedVelocityFromMatrix",
    "InfluenceMatrixTest.SupersonicMachClamped",
    "InfluenceMatrixTest.ComputeInducedVelocity",
    "InfluenceMatrixTest.ComputeSingleVortexVelocity",
    "InfluenceMatrixTest.YSymmetryInfluence",
    "InfluenceMatrixTest.DifferentComponentCoreRadius",
    "InfluenceMatrixTest.InducedVelocityWithYSymmetry",
    // mass_reader_test.cpp -> tests/mass_reader.rs
    "MassReaderTest.EmptyContent",
    "MassReaderTest.OnlyComments",
    "MassReaderTest.UnitDefinitions",
    "MassReaderTest.ConstantDefinitions",
    "MassReaderTest.SingleComponent",
    "MassReaderTest.ComponentWithProductsOfInertia",
    "MassReaderTest.MultipleComponents",
    "MassReaderTest.ScaleFactors",
    "MassReaderTest.OffsetFactors",
    "MassReaderTest.ScaleAndOffset",
    "MassReaderTest.TotalMass",
    "MassReaderTest.TotalMassWithUnits",
    "MassReaderTest.CenterOfGravity",
    "MassReaderTest.CenterOfGravityWeighted",
    "MassReaderTest.ComputeTotalProperties",
    "MassReaderTest.ParallelAxisTheorem",
    "MassReaderTest.EmptyComponentsGiveZeroMass",
    "MassReaderTest.CaseInsensitiveKeywords",
    "MassReaderTest.CompleteAVLStyleFile",
    // solver_test.cpp -> tests/solver_parity.rs
    "VortexLatticeSolverTest.SetAircraftGeneratesMesh",
    "VortexLatticeSolverTest.RebuildAICCreatesValidMatrix",
    "VortexLatticeSolverTest.SolveUnitResponsesComputes6Columns",
    "VortexLatticeSolverTest.SolveAtZeroAlphaGivesZeroCirculation",
    "VortexLatticeSolverTest.SolveAtPositiveAlphaGivesPositiveCirculation",
    "VortexLatticeSolverTest.CirculationLinearWithAlpha",
    "VortexLatticeSolverTest.FreestreamDirectionComputation",
    "VortexLatticeSolverTest.MachNumberChangeInvalidatesAIC",
    "VortexLatticeSolverTest.HighAspectRatioWingLiftSlope",
    "VortexLatticeSolverTest.EllipticLiftDistributionShape",
    "VortexLatticeSolverTest.YDuplicateProducesSymmetricCirculation",
    "VortexLatticeSolverTest.InducedVelocityAtControlPoints",
    "VortexLatticeSolverTest.GetMeshMutable",
    "VortexLatticeSolverTest.GetReferenceQuantities",
    "VortexLatticeSolverTest.SetReferencePoint",
    "VortexLatticeSolverTest.SetSymmetry",
    "VortexLatticeSolverTest.InitialStateInvalid",
    "VortexLatticeSolverTest.GetSolution",
    "VortexLatticeSolverTest.GetMatrices",
    "VortexLatticeSolverTest.FreestreamDirectionWithBeta",
    // vortex_lattice_test.cpp -> tests/vortex_lattice.rs
    "VortexLatticeTest.EmptyAircraft",
    "VortexLatticeTest.RectangularWingBasicCounts",
    "VortexLatticeTest.RectangularWingGeometry",
    "VortexLatticeTest.VortexQuarterChord",
    "VortexLatticeTest.ControlPointThreeQuarterChord",
    "VortexLatticeTest.NormalVectorsFlatWing",
    "VortexLatticeTest.TwistedWingNormals",
    "VortexLatticeTest.YDuplicateSymmetry",
    "VortexLatticeTest.TaperedWingChordDistribution",
    "VortexLatticeTest.SurfaceInfoConsistency",
    "VortexLatticeTest.StripVortexIndicesConsistent",
    "VortexLatticeTest.ClearAndRegenerate",
    "VortexLatticeTest.GetMeshInputNullByDefault",
    "VortexLatticeTest.MutableStripAccess",
    "VortexLatticeTest.GetVorticesMutable",
    "VortexLatticeTest.GetStripsMutable",
    "VortexLatticeTest.ControlSurfaceNamesEmpty",
    "VortexLatticeTest.ControlSurfaceInfluencesSize",
    "VortexLatticeTest.RecalculateNormals",
    "VortexLatticeTest.TwistOffsetEquivalentToGeometricTwist",
    // control_surface_test.cpp -> tests/control_surface.rs
    "ControlSurfaceDefinitionTest.NameDefaultsToEmpty",
    "ControlSurfaceDefinitionTest.SetValues",
    "SectionControlSurfaceTest.NoControlSurfacesByDefault",
    "SectionControlSurfaceTest.AddSingleControlSurface",
    "SectionControlSurfaceTest.AddMultipleControlSurfaces",
    "AircraftControlSurfaceTest.GetControlSurfaceNames",
    "AircraftControlSurfaceTest.GetNumControlSurfaces",
    "AircraftControlSurfaceTest.MultipleControlSurfaces",
    "MeshControlSurfaceTest.ControlSurfaceNamesInMesh",
    "MeshControlSurfaceTest.ControlInfluencesHaveCorrectSize",
    "MeshControlSurfaceTest.TrailingEdgePanelsHaveInfluence",
    "MeshControlSurfaceTest.LeadingEdgePanelsNoInfluence",
    "MeshControlSurfaceTest.InfluenceHingeAxisCorrect",
    "YDuplicateControlSurfaceTest.MirroredSurfaceCreated",
    "YDuplicateControlSurfaceTest.MirroredInfluencesExist",
    "YDuplicateControlSurfaceTest.MirroredControlSurfaceInfluences",
    "AvlControlSurfaceParsingTest.ParseBasicControlSurface",
    "AvlControlSurfaceParsingTest.ParseMultipleControlSurfaces",
    "AvlControlSurfaceParsingTest.ParseControlSurfaceHingeVector",
    "AvlControlSurfaceParsingTest.ParseControlSurfaceGain",
    "ControlSurfaceEffectivenessTest.ControlSurfacesPropagateToMesh",
    "ControlSurfaceEffectivenessTest.InfluenceMatrixBuiltWithControls",
    "ControlDerivativesTest.DefaultValues",
    "ControlDerivativesTest.SetValues",
    "FlightConditionControlTest.DefaultNoDeflections",
    "FlightConditionControlTest.SetControlDeflection",
    "StabilityDerivativesControlTest.ControlDerivativesVector",
    // eigenmode_test.cpp -> tests/eigenmode.rs
    "MassPropertiesTest.DefaultIsInvalid",
    "MassPropertiesTest.ValidWithPositiveValues",
    "MassPropertiesTest.InertiaTensorSymmetric",
    "MassPropertiesTest.InverseInertiaCorrect",
    "MassPropertiesTest.SetInertiaTensor",
    "FlightStateTest.BodyVelocityAtZeroAngles",
    "FlightStateTest.BodyVelocityWithAlpha",
    "FlightStateTest.DynamicPressure",
    "FlightStateTest.BodyVelocityWithBeta",
    "FlightStateTest.DefaultValues",
    "EigenmodeAnalysisTest.LongitudinalAnalysisProducesModes",
    "EigenmodeAnalysisTest.LateralDirectionalAnalysisProducesModes",
    "EigenmodeAnalysisTest.StableAircraftHasNegativeDamping",
    "EigenmodeAnalysisTest.ShortPeriodHigherFrequencyThanPhugoid",
    "EigenmodeAnalysisTest.EigenmodeGettersWork",
    "EigenmodeAnalysisTest.NonOscillatoryModeHasZeroPeriod",
    "EigenmodeAnalysisTest.DampingRatioComputation",
    "EigenmodeAnalysisTest.InvalidMassPropertiesFails",
    "EigenmodeAnalysisTest.InvalidDerivativesFails",
    "EigenmodeAnalysisTest.SystemMatrixHasCorrectSize",
    "EigenmodeAnalysisTest.FullAnalysisProducesMoreModes",
    "EigenmodeAnalysisTest.ModeTypeNames",
    "EigenmodeAnalysisTest.GetTimeConstant",
    "EigenmodeAnalysisTest.GetTimeConstantZeroDamping",
    "EigenmodeAnalysisTest.DampingRatioZeroFrequency",
    "EigenmodeAnalysisTest.UnstableModeDetected",
    "EigenmodeAnalysisTest.FindModesByTypeEmpty",
    "EigenmodeAnalysisTest.FindModesByTypeWithMatch",
    "EigenmodeAnalysisTest.GetDutchRollNotFound",
    "EigenmodeAnalysisTest.GetPhugoidNotFound",
    "EigenmodeAnalysisTest.GetShortPeriodNotFound",
    "EigenmodeAnalysisTest.IsStableAllStable",
    "EigenmodeAnalysisTest.IsStableOneUnstable",
    "EigenmodeAnalysisTest.IsStableEmpty",
    // geometry_test.cpp -> tests/geometry_parity.rs
    "AirfoilTest.FlatPlateCreation",
    "AirfoilTest.Naca0012Creation",
    "AirfoilTest.Naca2412Creation",
    "AirfoilTest.Naca4DigitFromString",
    "AirfoilTest.CamberSlopeAtLeadingEdge",
    "AirfoilExtraTest.InvalidNacaString",
    "AirfoilExtraTest.NacaCodeOutOfRange",
    "AirfoilExtraTest.NacaZeroThickness",
    "AirfoilExtraTest.EvaluateThicknessNaca0012",
    "AirfoilExtraTest.EvaluateCamberClamping",
    "AirfoilExtraTest.InvalidAirfoilReturnsZero",
    "AirfoilExtraTest.SetName",
    "AirfoilExtraTest.CreateFromCoordinatesCamberLine",
    "AirfoilExtraTest.CreateFromCoordinatesMismatch",
    "AirfoilExtraTest.CreateFromCoordinatesTooFew",
    "AirfoilExtraTest.CreateFromCoordinatesFullAirfoil",
    "AirfoilExtraTest.CreateFromCoordinateFileNotFound",
    "AirfoilExtraTest.CamberXY",
    "SectionTest.DefaultConstruction",
    "SectionTest.ParameterizedConstruction",
    "SectionTest.AirfoilAssignment",
    "SectionTest.ControlSurface",
    "SurfaceTest.BasicConstruction",
    "SurfaceTest.AddSections",
    "SurfaceTest.PlanformAreaRectangular",
    "SurfaceTest.PlanformAreaTapered",
    "SurfaceTest.YDuplicateDoublesArea",
    "SurfaceExtraTest.MeanAerodynamicChordRectangular",
    "SurfaceExtraTest.MeanAerodynamicChordTapered",
    "SurfaceExtraTest.MeanAerodynamicChordSingleSection",
    "SurfaceExtraTest.MeanAerodynamicChordEmpty",
    "SurfaceExtraTest.ComputeSpan",
    "SurfaceExtraTest.ComputeSpanYDuplicate",
    "SurfaceExtraTest.ComputeSpanEmpty",
    "SurfaceExtraTest.PlanformAreaSingleSection",
    "SurfaceExtraTest.HasViscousDragDefault",
    "SurfaceExtraTest.HasViscousDragNonZero",
    "SurfaceExtraTest.ChordwiseSpacingNumeric",
    "SurfaceExtraTest.SpanwiseSpacingNumeric",
    "SurfaceExtraTest.ClearYDuplicate",
    "SurfaceExtraTest.AerodynamicOptions",
    "SurfaceExtraTest.CoordinateTransforms",
    "SurfaceExtraTest.ComponentIndex",
    "SurfaceExtraTest.SectionAccessors",
    "AircraftTest.DefaultValues",
    "AircraftTest.AddSurface",
    "AircraftTest.FindSurfaceByName",
    "AircraftTest.ReferenceQuantities",
    "AircraftTest.FlowConditions",
    "AircraftTest.Symmetry",
    "AircraftTest.Bodies",
    "AircraftTest.SurfaceAccessors",
    "AircraftTest.MoveSurface",
    "AircraftTest.ControlSurfaceNames",
    "AircraftTest.DesignVariableNames",
    "AircraftTest.ComputeTotalVortexCount",
    "AircraftTest.ComputeTotalVortexCountYDuplicate",
    "AircraftTest.ComputeTotalStripCount",
    "AircraftTest.ComputeTotalStripCountYDuplicate",
    "AircraftTest.ComputeTotalPlanformArea",
    "AircraftTest.ComputeTotalSpan",
    "AvlReaderTest.ParseSimpleGeometry",
    "AvlReaderTest.ParseWithNacaAirfoil",
    "AvlReaderTest.ParseWithYDuplicate",
    "AvlReaderTest.ParseWithControlSurface",
    "AvlReaderExtraTest.ReadFileNotFound",
    "AvlReaderExtraTest.ParseEmptyContent",
    "AvlReaderExtraTest.ParseWithComments",
    "AvlReaderExtraTest.ParseWithProfileDrag",
    "AvlReaderExtraTest.ParseWithSymmetry",
    "AvlReaderExtraTest.ParseSurfaceModifiers",
    "AvlReaderExtraTest.ParseSectionSpanwise",
    "AvlReaderExtraTest.ParseWithDesignVariable",
    "AvlReaderExtraTest.ParseWithBody",
    "AvlReaderExtraTest.ParseBodyAtEnd",
    "AvlReaderExtraTest.ParseVerbose",
    "AvlReaderExtraTest.ParseMomentReference",
    "AvlReaderExtraTest.ParseMultipleSurfaces",
    // force_calculator_test.cpp -> tests/force_calculator.rs
    "ForceCalculatorTest.ZeroAlphaGivesZeroForces",
    "ForceCalculatorTest.PositiveAlphaGivesPositiveLift",
    "ForceCalculatorTest.NegativeAlphaGivesNegativeLift",
    "ForceCalculatorTest.LiftCoefficientMatchesCirculation",
    "ForceCalculatorTest.SpanEfficiencyReasonable",
    "ForceCalculatorTest.InducedDragQuadraticWithLift",
    "ForceCalculatorTest.SurfaceLoadsAddUpToTotal",
    "ForceCalculatorTest.StripLoadsHaveCorrectCount",
    "ForceCalculatorTest.PitchingMomentSign",
    "ForceCalculatorTest.AspectRatioCalculation",
    "ForceCalculatorTest.IdealInducedDragComputation",
    // stability_derivatives_test.cpp -> tests/force_calculator.rs
    "StabilityDerivativesTest.LiftCurveSlopePositive",
    "StabilityDerivativesTest.LiftCurveSlopeReasonable",
    "StabilityDerivativesTest.ZeroAlphaSymmetry",
    "StabilityDerivativesTest.RollDampingNegative",
    "StabilityDerivativesTest.PitchDampingFinite",
    "StabilityDerivativesTest.YawDampingNonPositive",
    "StabilityDerivativesTest.BodyAxisDerivativesValid",
    "StabilityDerivativesTest.DerivativesAtDifferentAlpha",
    "StabilityDerivativesTest.StaticMarginComputed",
    "StabilityDerivativesTest.ReferenceQuantitiesPreserved",
    "StabilityDerivativesTest.HighAspectRatioRollDamping",
    // section_test.cpp -> tests/section.rs
    "SectionUnitTest.DefaultConstructor",
    "SectionUnitTest.ParameterizedConstructor",
    "SectionUnitTest.ConstructorWithAirfoil",
    "SectionUnitTest.GetAirfoilWithoutCustomReturnsFlat",
    "SectionUnitTest.GetAirfoilMutableCreatesFlat",
    "SectionUnitTest.GetAirfoilWithCustomReturnsCustom",
    "SectionUnitTest.SetAirfoil",
    "SectionUnitTest.SetLeadingEdgePosition",
    "SectionUnitTest.SetChord",
    "SectionUnitTest.SetTwistAngle",
    "SectionUnitTest.TwistAngleRadians",
    "SectionUnitTest.ClAlphaFactor",
    "SectionUnitTest.TwistOffsetDefaultsToZero",
    "SectionUnitTest.TwistOffsetAddsToGeometricTwist",
    "SectionUnitTest.AddControlSurface",
    "SectionUnitTest.AddDesignVariable",
    "SectionUnitTest.SpanwisePanelCount",
    "SectionUnitTest.SpanwiseSpacingEnum",
    "SectionUnitTest.SpanwiseSpacingNumericUniform",
    "SectionUnitTest.SpanwiseSpacingNumericCosine",
    "SectionUnitTest.SpanwiseSpacingNumericSine",
    "SectionUnitTest.ViscousPolar",
    // flight_constraints_test.cpp -> tests/flight_constraints.rs
    "OperatingVariablesTest.DefaultNumVariables",
    "OperatingVariablesTest.NumVariablesWithControls",
    "OperatingVariablesTest.DefaultValues",
    "TrimConstraintsTest.DefaultAllFree",
    "TrimConstraintsTest.ConstraintCountWithControls",
    "TrimConstraintsTest.AllConstraintTypesCount",
    "TrimPresetsTest.LevelFlightRatesFixed",
    "TrimPresetsTest.LevelTrimmedFlight",
    "TrimPresetsTest.SteadyTurn",
    "TrimPresetsTest.SteadySideslipAllFields",
    "VariableConstraintTest.DefaultState",
];

/// Frozen cases that cannot be ported because the feature they cover no longer
/// exists, deliberately.
///
/// These are not a porting backlog. `ApparentMass` was carried by v0.2.0 but
/// never read by `eigenmode_analysis.cpp`, so there was no behaviour to
/// reproduce — the struct and its setter have been removed rather than kept as
/// dead API. The two cases matched that dead surface: `DefaultState` asserted
/// the zero defaults of a struct that fed nothing, and `SetApparentMass`
/// asserted only, in its own words, that the setter "doesn't crash".
///
/// Membership is verified against the manifest, and a name may not appear in
/// both this list and [`PORTED`].
const UNPORTABLE: &[&str] = &[
    "ApparentMassTest.DefaultState",
    "EigenmodeAnalysisTest.SetApparentMass",
    // vortex_lattice_test.cpp: `VortexLattice::{Get,Set}VortexCoreRadius*` was a
    // second pair of core-radius knobs never passed to assembly — setting them
    // left the AIC bit-identical — and has been removed. `MutableVortexAccess`
    // asserted in-place panel mutation through a public reference; `vortices_mut`
    // is crate-private so the lattice cannot be desynchronised from its caches.
    "VortexLatticeTest.CoreRadiusFractions",
    "VortexLatticeTest.MutableVortexAccess",
];

const CATEGORIES: &[&str] = &[
    "aero",
    "dynamics",
    "geometry",
    "gradients",
    "io",
    "solver",
    "trim",
    "visualization",
];

fn manifest() -> Manifest {
    serde_json::from_str(include_str!("fixtures/parity/ctest-manifest.json"))
        .expect("frozen parity manifest must parse")
}

#[test]
fn manifest_matches_its_recorded_provenance() {
    let manifest = manifest();
    assert_eq!(manifest.provenance.tag, "v0.2.0");
    assert_eq!(manifest.provenance.commit, "556b0ea");
    assert_eq!(manifest.cases.len(), manifest.provenance.rust_parity_cases);
    assert_eq!(
        manifest.provenance.ctest_cases,
        manifest.provenance.rust_parity_cases + manifest.provenance.deferred_philote_cases.len(),
        "recorded C++ case count must equal ported-scope cases plus deferred Philote cases"
    );
}

#[test]
fn manifest_case_inventory_is_well_formed() {
    let manifest = manifest();
    let mut seen = HashSet::new();
    for case in &manifest.cases {
        assert!(
            seen.insert(case.cpp_test.clone()),
            "duplicate case in manifest: {}",
            case.cpp_test
        );
        assert!(
            case.cpp_test.contains('.'),
            "case name must be `Suite.Case`: {}",
            case.cpp_test
        );
        assert!(
            CATEGORIES.contains(&case.category.as_str()),
            "unknown category {} on {}",
            case.category,
            case.cpp_test
        );
        assert!(
            case.tolerance >= 0.0 && case.tolerance.is_finite(),
            "tolerance must be finite and non-negative on {}",
            case.cpp_test
        );
        assert!(
            case.enabled,
            "disabled case left in manifest: {}",
            case.cpp_test
        );
    }
}

#[test]
fn ported_cases_exist_in_the_manifest() {
    let manifest = manifest();
    let known: HashSet<&str> = manifest.cases.iter().map(|c| c.cpp_test.as_str()).collect();
    for name in PORTED.iter().chain(UNPORTABLE) {
        assert!(
            known.contains(name),
            "PORTED/UNPORTABLE names a case that is not in the frozen manifest: {name}"
        );
    }
    let ported: HashSet<&str> = PORTED.iter().copied().collect();
    for name in UNPORTABLE {
        assert!(
            !ported.contains(name),
            "{name} is listed as both ported and unportable"
        );
    }
}

/// Guards the porting gap so it cannot be misreported in either direction.
///
/// This test failing means the number of ported or unportable C++ cases
/// changed. That is expected when work lands — update the counts deliberately
/// and say so in `CHANGELOG.md`.
#[test]
fn ported_case_count_is_recorded() {
    let manifest = manifest();
    let ported = PORTED.len();
    let unportable = UNPORTABLE.len();
    let total = manifest.cases.len();
    assert_eq!(
        ported, 353,
        "{ported} of {total} frozen C++ cases are individually ported; \
         update this count when that changes"
    );
    assert_eq!(
        unportable, 4,
        "{unportable} of {total} frozen C++ cases cover removed behaviour; \
         update this count when that changes"
    );
    assert_eq!(
        total - ported - unportable,
        131,
        "the remaining porting backlog changed; update this count deliberately"
    );
}
