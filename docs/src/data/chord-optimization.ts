// Figures 2 and 3 are drawn from one run of `examples/chord_optimization.py`,
// the OpenMDAO driver, whose own output is committed at
// `examples/chord_optimization_results.json`. These are that run's numbers, so
// the figures are the result as solved rather than an illustration of it.
//
// Regenerate by running the example and re-exporting this file from its JSON.

export const problem = {
    span: 10.0,
    area: 10.0,
    aspectRatio: 10.0,
    targetLift: 0.5,
    stations: 21,
    controlPoints: 6,
    splineOrder: 4,
    chordBounds: [0.05, 2.5] as const,
    spanwisePanels: 20,
    chordwisePanels: 4,
};

/** Semispan station positions the spline is evaluated at. */
export const stations = [0.000000, 0.250000, 0.500000, 0.750000, 1.000000, 1.250000, 1.500000, 1.750000, 2.000000, 2.250000, 2.500000, 2.750000, 3.000000, 3.250000, 3.500000, 3.750000, 4.000000, 4.250000, 4.500000, 4.750000, 5.000000];

/** Chord at each station: the rectangular start, the optimum, and the ellipse. */
export const initialChords = [1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000];
export const optimizedChords = [1.239215, 1.237756, 1.233999, 1.227639, 1.218369, 1.205886, 1.189882, 1.170055, 1.146193, 1.118231, 1.086116, 1.049793, 1.009209, 0.964312, 0.914851, 0.858298, 0.790653, 0.707893, 0.605992, 0.480927, 0.328674];
export const ellipticalChords = [1.273240, 1.271647, 1.266857, 1.258834, 1.247515, 1.232809, 1.214593, 1.192707, 1.166943, 1.137039, 1.102658, 1.063365, 1.018592, 0.967578, 0.909275, 0.842169, 0.763944, 0.670721, 0.554992, 0.397569, 0.000000];

/** IPOPT primal infeasibility (inf_pr) at each iteration: feasibility. */
export const feasibility = [2.220e-16, 4.880e-05, 7.840e-03, 2.640e-03, 7.240e-04, 2.290e-03, 6.320e-04, 3.510e-05, 1.690e-05, 1.390e-05, 6.180e-06, 9.860e-06, 4.010e-06, 6.080e-07, 1.320e-07, 6.040e-09, 4.110e-09, 4.980e-08, 1.980e-07, 5.870e-08, 8.790e-06, 5.670e-06, 5.000e-07, 3.000e-07, 2.510e-10, 2.220e-09, 5.150e-09, 3.620e-09, 1.440e-12, 1.680e-13];

/** IPOPT dual infeasibility (inf_du) at each iteration: optimality. */
export const optimality = [5.060e-02, 1.680e-02, 9.520e-03, 1.720e-02, 9.480e-03, 9.110e-03, 4.310e-03, 1.340e-03, 6.640e-04, 4.580e-04, 5.330e-04, 7.740e-04, 1.500e-04, 4.790e-05, 3.190e-05, 1.950e-05, 3.870e-05, 2.900e-05, 8.990e-05, 6.590e-05, 4.380e-04, 5.790e-05, 1.510e-04, 1.050e-05, 6.410e-06, 2.650e-06, 1.840e-05, 5.730e-07, 1.200e-07, 6.610e-09];
