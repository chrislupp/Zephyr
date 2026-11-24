import type { ReactNode } from "react";
import Head from "@docusaurus/Head";
import Link from "@docusaurus/Link";
import CodeBlock from "@theme/CodeBlock";
import Layout from "@theme/Layout";
import ContentsRail from "@site/src/components/ContentsRail";
import LatticeFigure from "@site/src/components/LatticeFigure";
import { docs, sheetNumber } from "@site/src/data/contents";
import styles from "./index.module.css";

// The contents page. It carries no sheet number of its own; its rail lists
// both books, numbered the way each book's own rail numbers them.

type RouteSheet = { label: string; to: string };

const routes: { title: string; body: string; sheets: RouteSheet[] }[] = [
    {
        title: "Run a first analysis",
        body: "Add the crate, describe a wing, and read its loads.",
        sheets: [
            { label: "Installation", to: "/docs/getting-started/installation" },
            { label: "Quick Start", to: "/docs/getting-started/quickstart" },
            { label: "Conventions", to: "/docs/getting-started/conventions" },
            { label: "Basic Analysis", to: "/docs/examples/basic-analysis" },
        ],
    },
    {
        title: "Wire it into an optimizer",
        body: "Analytical and adjoint gradients, and Philote servers for MDO frameworks.",
        sheets: [
            { label: "Analytical and Adjoint Gradients",
                to: "/docs/guide/automatic-differentiation",
            },
            { label: "Selective Analysis", to: "/docs/guide/selective-analysis" },
            { label: "Coupled Interfaces", to: "/docs/guide/coupled-interfaces" },
            { label: "MDO Integration", to: "/docs/guide/mdo-integration" },
        ],
    },
];

const conventions: [string, string][] = [
    ["Axes", "Right-handed body frame: x aft, y starboard, z up."],
    ["Degrees", "Section twist and surface incidence."],
    [
        "Radians",
        "Flight angles, angular rates, control deflections, trim variables, and derivatives.",
    ],
    [
        "Units",
        "Length, area, mass, density, and dynamic pressure in one consistent system chosen by the caller.",
    ],
    ["Precision", "Double precision throughout."],
];

export default function Home(): ReactNode {
    return (
        <Layout
            title="Documentation"
            description="Zephyr is a vortex-lattice-method library in safe Rust: aircraft geometry, aerodynamic loads, stability derivatives, trim, eigenmodes, and analytical and adjoint gradients."
        >
            <Head>
                <html className="zm-has-rail" />
            </Head>
            <div className={styles.sheet}>
                <ContentsRail withApi />
                <main className={styles.main}>
                    <div className={styles.intro}>
                        <h1 className={styles.name}>Zephyr</h1>
                        <p className={styles.lede}>
                            A vortex-lattice-method library in safe Rust for aircraft
                            geometry, aerodynamic loads, stability and control
                            derivatives, trim, eigenmodes, and analytical and adjoint
                            gradients. It reads AVL geometry and mass files and writes
                            VTK.
                        </p>
                        <CodeBlock language="bash">cargo add zephyr-vlm</CodeBlock>
                    </div>
                    <ol className={styles.routes}>
                        {routes.map((route) => (
                            <li key={route.title} className={styles.route}>
                                <h2>{route.title}</h2>
                                <p>{route.body}</p>
                                <ul>
                                    {route.sheets.map((sheet) => (
                                        <li key={sheet.to}>
                                            <Link to={sheet.to}>
                                                <span className={styles.routeNumber}>
                                                    {sheetNumber(docs, sheet.to)}
                                                </span>
                                                {sheet.label}
                                            </Link>
                                        </li>
                                    ))}
                                </ul>
                            </li>
                        ))}
                    </ol>

                    <figure className={styles.figure}>
                        <div className={styles.figureScroll}>
                            <LatticeFigure />
                        </div>
                        <figcaption>
                            <b>Figure 1.</b> The Quick Start wing as Zephyr panels it: 4
                            chordwise by 16 spanwise panels, cosine spaced. Bound
                            vortices in red, control points as dots.
                            <span className={styles.narrowNote}>
                                {" "}
                                The starboard half is shown.
                            </span>
                        </figcaption>
                    </figure>

                    <section className={styles.conventions}>
                        <h2>Conventions</h2>
                        <dl>
                            {conventions.map(([term, definition]) => (
                                <div key={term}>
                                    <dt>{term}</dt>
                                    <dd>{definition}</dd>
                                </div>
                            ))}
                        </dl>
                    </section>
                </main>
            </div>
        </Layout>
    );
}
