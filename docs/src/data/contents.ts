// The contents of the manual, as the index sheet prints them. Keep this in
// step with `sidebars.ts`: same sections, same order. Sheet numbers are not
// stored here; the page counts them the same way the rail does. Every `to` is
// checked by the build (`onBrokenLinks: "throw"`).

export type Sheet = { label: string; to: string };
export type Section = { label: string; to?: string; sheets: Sheet[] };

export const docs: Section[] = [
    {
        label: "Getting Started",
        sheets: [
            { label: "Installation", to: "/docs/getting-started/installation" },
            { label: "Quick Start", to: "/docs/getting-started/quickstart" },
            { label: "Conventions", to: "/docs/getting-started/conventions" },
        ],
    },
    {
        label: "Guide",
        sheets: [
            { label: "VLM Theory", to: "/docs/guide/vlm-theory" },
            { label: "Geometry", to: "/docs/guide/geometry" },
            { label: "Control Surfaces", to: "/docs/guide/control-surfaces" },
            { label: "Solver Workflow", to: "/docs/guide/solver-workflow" },
            { label: "Loads and Coefficients", to: "/docs/guide/loads" },
            {
                label: "Stability and Control Derivatives",
                to: "/docs/guide/stability-analysis",
            },
            { label: "Trim", to: "/docs/guide/trim" },
            { label: "Dynamic Modes", to: "/docs/guide/dynamic-modes" },
            { label: "AVL Input Files", to: "/docs/guide/avl-files" },
            {
                label: "Analytical and Adjoint Gradients",
                to: "/docs/guide/automatic-differentiation",
            },
            {
                label: "Selective Analysis",
                to: "/docs/guide/selective-analysis",
            },
            {
                label: "Coupled Interfaces",
                to: "/docs/guide/coupled-interfaces",
            },
            { label: "Performance", to: "/docs/guide/performance" },
            { label: "Validity and Limits", to: "/docs/guide/limitations" },
            { label: "MDO Integration", to: "/docs/guide/mdo-integration" },
        ],
    },
    {
        label: "Examples",
        sheets: [
            { label: "Basic Analysis", to: "/docs/examples/basic-analysis" },
            { label: "Trim Flight", to: "/docs/examples/trim-flight" },
            { label: "Stability Modes", to: "/docs/examples/stability-modes" },
            { label: "VTK Visualization", to: "/docs/examples/vtk-visualization" },
            { label: "Control Surfaces", to: "/docs/examples/control-surfaces" },
            { label: "Mesh Trim", to: "/docs/examples/mesh-trim" },
            {
                label: "Selective Analysis",
                to: "/docs/examples/selective-analysis",
            },
            {
                label: "Philote Optimization",
                to: "/docs/examples/philote-optimization",
            },
        ],
    },
    { label: "Copyright and License", to: "/public-release", sheets: [] },
];

export const api: Section[] = [
    { label: "API Reference", to: "/docs/api", sheets: [] },
    {
        label: "Geometry",
        sheets: [
            { label: "Aircraft", to: "/docs/api/aircraft" },
            { label: "Surface", to: "/docs/api/surface" },
            { label: "Section", to: "/docs/api/section" },
            { label: "MeshInput", to: "/docs/api/mesh-input" },
        ],
    },
    {
        label: "Solver",
        sheets: [
            { label: "VortexLatticeSolver", to: "/docs/api/vortex-lattice-solver" },
            { label: "ForceCalculator", to: "/docs/api/force-calculator" },
            { label: "Analysis", to: "/docs/api/analysis" },
        ],
    },
    {
        label: "Analysis",
        sheets: [
            {
                label: "StabilityDerivativeCalculator",
                to: "/docs/api/stability-derivatives",
            },
            { label: "TrimSolver", to: "/docs/api/trim-solver" },
            { label: "MeshTrimSolver", to: "/docs/api/mesh-trim-solver" },
            { label: "EigenmodeAnalysis", to: "/docs/api/eigenmode-analysis" },
        ],
    },
    {
        label: "Gradients",
        sheets: [
            { label: "GradientCalculator", to: "/docs/api/gradient-calculator" },
            { label: "StripLoadAdjoint", to: "/docs/api/strip-load-adjoint" },
        ],
    },
    {
        label: "I/O",
        sheets: [
            { label: "AVL and Mass Readers", to: "/docs/api/avl-reader" },
            { label: "VtkExporter", to: "/docs/api/visualization" },
        ],
    },
    {
        label: "Errors",
        sheets: [{ label: "ZephyrError", to: "/docs/api/errors" }],
    },
];

/** The sheet number of `to` within a book: "4" or "2.1"; empty if absent. */
export function sheetNumber(sections: Section[], to: string, prefix = ""): string {
    for (const [i, section] of sections.entries()) {
        if (section.to === to) {
            return `${prefix}${i + 1}`;
        }
        const j = section.sheets.findIndex((sheet) => sheet.to === to);
        if (j >= 0) {
            return `${prefix}${i + 1}.${j + 1}`;
        }
    }
    return "";
}
