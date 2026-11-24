import type { SidebarsConfig } from "@docusaurus/plugin-content-docs";

const sidebars: SidebarsConfig = {
    docsSidebar: [
        {
            type: "category",
            label: "Getting Started",
            items: [
                "getting-started/installation",
                "getting-started/quickstart",
                "getting-started/conventions",
            ],
        },
        {
            type: "category",
            label: "Guide",
            items: [
                "guide/vlm-theory",
                "guide/geometry",
                "guide/control-surfaces",
                "guide/solver-workflow",
                "guide/loads",
                "guide/stability-analysis",
                "guide/trim",
                "guide/dynamic-modes",
                "guide/avl-files",
                "guide/automatic-differentiation",
                "guide/selective-analysis",
                "guide/coupled-interfaces",
                "guide/performance",
                "guide/limitations",
                "guide/mdo-integration",
            ],
        },
        {
            type: "category",
            label: "Examples",
            items: [
                "examples/basic-analysis",
                "examples/trim-flight",
                "examples/stability-modes",
                "examples/vtk-visualization",
                "examples/control-surfaces",
                "examples/mesh-trim",
                "examples/selective-analysis",
                "examples/chord-optimization",
                "examples/philote-optimization",
            ],
        },
        {
            type: "link",
            label: "Copyright and License",
            href: "/public-release",
        },
    ],
    apiSidebar: [
        "api/index",
        {
            type: "category",
            label: "Geometry",
            items: [
                "api/aircraft",
                "api/surface",
                "api/section",
                "api/mesh-input",
            ],
        },
        {
            type: "category",
            label: "Solver",
            items: [
                "api/vortex-lattice-solver",
                "api/force-calculator",
                "api/analysis",
            ],
        },
        {
            type: "category",
            label: "Analysis",
            items: [
                "api/stability-derivatives",
                "api/trim-solver",
                "api/mesh-trim-solver",
                "api/eigenmode-analysis",
            ],
        },
        {
            type: "category",
            label: "Gradients",
            items: [
                "api/gradient-calculator",
                "api/strip-load-adjoint",
            ],
        },
        {
            type: "category",
            label: "I/O",
            items: [
                "api/avl-reader",
                "api/visualization",
            ],
        },
        {
            type: "category",
            label: "Errors",
            items: [
                "api/errors",
            ],
        },
    ],
};

export default sidebars;
