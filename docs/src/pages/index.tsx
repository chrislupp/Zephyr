import { useEffect, useRef, useState, type ReactNode } from "react";
import Head from "@docusaurus/Head";
import Link from "@docusaurus/Link";
import useBaseUrl from "@docusaurus/useBaseUrl";
import CodeBlock from "@theme/CodeBlock";
import Layout from "@theme/Layout";
import {
    ellipticalChords,
    initialChords,
    optimizedChords,
    stations,
} from "@site/src/data/chord-optimization";
import styles from "./index.module.css";

// The start page: a post-processor's result window. Zephyr's own solution on
// the uCRM transport wing fills the first screen; below it, three things the
// library does, each shown by its own output, and one short listing. The
// renders come from `docs/scripts/render_loading.py`; the numbers are the
// solver's for that run.

const INSTALL = "cargo add zephyr-vlm";

/** The uCRM run the hero shows: tests/fixtures/avl_reference/ucrm_17_40.avl. */
const HERO = { alpha: "2.5°", panels: 680, lift: "0.1723", drag: "0.00340" };

function InstallButton({ className }: { className: string }): ReactNode {
    const [copied, setCopied] = useState(false);
    const reset = useRef<ReturnType<typeof setTimeout> | undefined>(undefined);
    useEffect(() => () => clearTimeout(reset.current), []);
    const copy = async (): Promise<void> => {
        try {
            await navigator.clipboard.writeText(INSTALL);
            setCopied(true);
            clearTimeout(reset.current);
            reset.current = setTimeout(() => setCopied(false), 1600);
        } catch {
            // No clipboard access: the command stays on screen to select.
        }
    };
    return (
        <button type="button" className={className} onClick={copy}>
            <code>{INSTALL}</code>
            <span className={styles.copyState} aria-live="polite">
                {copied ? "Copied" : "Copy"}
            </span>
        </button>
    );
}

function Hero(): ReactNode {
    const mesh = useBaseUrl("/img/renders/ucrm-mesh.svg");
    const loading = useBaseUrl("/img/renders/ucrm-loading.svg");
    return (
        <header className={styles.hero}>
            <figure className={styles.render}>
                <img src={mesh} alt="" className={styles.mesh} width={1400} height={720} />
                <img
                    src={loading}
                    alt={`The uCRM transport wing as Zephyr panels and solves it: ${HERO.panels} panels at an angle of attack of ${HERO.alpha}, coloured by panel loading, with the trailing vortex legs running aft.`}
                    className={styles.solution}
                    width={1400}
                    height={720}
                />
            </figure>
            <p className={styles.caption}>
                uCRM wing, {HERO.panels} panels
                <br />
                α = {HERO.alpha}, solved by Zephyr
            </p>
            <div className={styles.copy}>
                <h1 className={styles.title}>Vortex-lattice aerodynamics, with exact gradients.</h1>
                <p className={styles.lede}>
                    Loads, stability derivatives, trim and eigenmodes in safe Rust,
                    differentiated by hand for design optimization.
                </p>
                <div className={styles.actions}>
                    <Link className={styles.primary} to="/docs/getting-started/quickstart">
                        Get started
                    </Link>
                    <InstallButton className={styles.install} />
                </div>
            </div>
            <dl className={styles.readouts}>
                <div>
                    <dt>
                        C<sub>L</sub>
                    </dt>
                    <dd>{HERO.lift}</dd>
                </div>
                <div>
                    <dt>
                        C<sub>Di</sub>
                    </dt>
                    <dd>{HERO.drag}</dd>
                </div>
                <div>
                    <dt>Checked against</dt>
                    <dd>
                        <Link to="/docs/guide/limitations#what-has-been-checked">OptVL 2.1.1</Link>
                    </dd>
                </div>
                <div className={styles.scale}>
                    <dt>
                        Panel loading ΔC<sub>p</sub>
                    </dt>
                    <dd>
                        <i className={styles.ramp} aria-hidden="true" />
                        <span className={styles.ticks}>
                            <span>0</span>
                            <span>0.5</span>
                            <span>≥ 1.0</span>
                        </span>
                    </dd>
                </div>
            </dl>
        </header>
    );
}

// ── Capability visuals ──────────────────────────────────────────────────────

/** The chord-optimization planforms, starboard half, at one scale. */
function Planform(): ReactNode {
    const w = 360;
    const k = w / 5;
    const top = 6;
    const root = Math.max(...ellipticalChords);
    // Quarter-chord line held straight, as the example holds it.
    const le = (c: number): number => top + 0.25 * (root - c) * k;
    const outline = (chords: number[]): string => {
        const lead = chords.map((c, j) => `${(stations[j] * k).toFixed(1)},${le(c).toFixed(1)}`);
        const trail = chords
            .map((c, j) => `${(stations[j] * k).toFixed(1)},${(le(c) + c * k).toFixed(1)}`)
            .reverse();
        return `M${[...lead, ...trail].join(" L")}`;
    };
    return (
        <svg viewBox={`-4 0 ${w + 8} ${top * 2 + root * k}`} className={styles.planform} aria-hidden="true">
            <path d={outline(initialChords)} className={styles.pfStart} />
            <path d={outline(ellipticalChords)} className={styles.pfTarget} />
            <path d={outline(optimizedChords)} className={styles.pfSolved} />
        </svg>
    );
}

const fixtures: [string, number][] = [
    ["rectangular", 8],
    ["rectangular_claf", 8],
    ["rectangular_sine", 8],
    ["wing_tail", 44],
    ["ucrm_5_8", 40],
    ["ucrm_17_40", 680],
];

type Capability = {
    title: string;
    body: string;
    links: { label: string; to: string }[];
    visual: ReactNode;
    note: string;
};

function Capabilities(): ReactNode {
    const wingTail = useBaseUrl("/img/renders/wingtail-loading.svg");
    const items: Capability[] = [
        {
            title: "Analysis",
            body: "Near- and far-field loads, stability and control derivatives, trim, and dynamic modes, from AVL geometry or a mesh you supply.",
            links: [
                { label: "Loads", to: "/docs/guide/loads" },
                { label: "Stability", to: "/docs/guide/stability-analysis" },
                { label: "Trim", to: "/docs/guide/trim" },
            ],
            visual: (
                <img
                    src={wingTail}
                    alt="A wing and tail solved by Zephyr, 44 panels at 5 degrees angle of attack, coloured by panel loading."
                    className={styles.wingTail}
                    width={800}
                    height={560}
                    loading="lazy"
                />
            ),
            note: "Wing and tail, 44 panels, α = 5°, M = 0.3",
        },
        {
            title: "Gradients for MDO",
            body: "Every derivative is analytical. An adjoint costs one transposed solve, and Philote serves the mesh and the aerodynamics to OpenMDAO.",
            links: [
                { label: "Gradients", to: "/docs/guide/automatic-differentiation" },
                { label: "MDO integration", to: "/docs/guide/mdo-integration" },
            ],
            visual: <Planform />,
            note: "Chord optimization: start, IPOPT result, ellipse",
        },
        {
            title: "Validated",
            body: "Checked against independent OptVL 2.1.1 fixtures: coefficients to 2e-8 and derivatives to 2e-7, up to the 680-panel uCRM.",
            links: [{ label: "Validity and limits", to: "/docs/guide/limitations" }],
            visual: (
                <ul className={styles.fixtures} aria-label="OptVL validation cases">
                    {fixtures.map(([name, panels]) => (
                        <li key={name}>
                            <code>{name}</code>
                            <span>{panels} panels</span>
                        </li>
                    ))}
                </ul>
            ),
            note: "tests/fixtures/avl_reference/",
        },
    ];
    return (
        <section className={styles.section} aria-labelledby="capabilities">
            <h2 id="capabilities" className={styles.sectionTitle}>
                The AVL-style workflow, and the optimizer around it.
            </h2>
            <div className={styles.capabilities}>
                {items.map((item) => (
                    <article key={item.title} className={styles.capability}>
                        <div className={styles.visual}>{item.visual}</div>
                        <p className={styles.note}>{item.note}</p>
                        <h3>{item.title}</h3>
                        <p>{item.body}</p>
                        <p className={styles.links}>
                            {item.links.map((link) => (
                                <Link key={link.to} to={link.to}>
                                    {link.label}
                                </Link>
                            ))}
                        </p>
                    </article>
                ))}
            </div>
        </section>
    );
}

const listing = `let mut wing = Surface::new("Wing");
wing.set_chordwise_panel_count(4);
wing.set_spanwise_panel_count(16);
wing.add_section(Section::new(Vec3::new(0.0, -5.0, 0.0), 1.0, 0.0));
wing.add_section(Section::new(Vec3::new(0.0, 5.0, 0.0), 1.0, 0.0));

let mut aircraft = Aircraft::new("Rectangular wing");
aircraft.set_reference_area(10.0);
aircraft.set_reference_chord(1.0);
aircraft.set_reference_span(10.0);
aircraft.add_surface(wing);

let condition = FlightCondition {
    alpha: 5.0 * DEGREES_TO_RADIANS, ..Default::default()
};

let mut solver = VortexLatticeSolver::new();
solver.set_aircraft(&aircraft)?;
solver.solve(&condition)?;

let lift = TotalDerivative::compute(&solver, &condition, AeroOutput::Lift)?;
println!("dCL/dalpha = {}", lift.by_alpha);`;

function FirstAnalysis(): ReactNode {
    return (
        <section className={`${styles.section} ${styles.code}`} aria-labelledby="first-analysis">
            <div className={styles.codeText}>
                <h2 id="first-analysis" className={styles.sectionTitle}>
                    A wing, a solve, and its derivative.
                </h2>
                <p>
                    Describe the geometry by sections or read an AVL file, solve, and ask
                    for loads or gradients from the same factorization. The Quick Start
                    walks through it.
                </p>
                <div className={styles.actions}>
                    <Link className={styles.primaryAlt} to="/docs/getting-started/quickstart">
                        Read the Quick Start
                    </Link>
                    <Link className={styles.secondary} to="/docs/api">
                        API reference
                    </Link>
                </div>
            </div>
            <div className={styles.codeBlock}>
                <CodeBlock language="rust">{listing}</CodeBlock>
            </div>
        </section>
    );
}

export default function Home(): ReactNode {
    return (
        <Layout
            title="Vortex-lattice aerodynamics with exact gradients"
            description="Zephyr is a vortex-lattice-method library in safe Rust: aircraft geometry, aerodynamic loads, stability derivatives, trim, eigenmodes, and analytical and adjoint gradients for design optimization."
        >
            <Head>
                <html className="z-home" />
            </Head>
            <main>
                <Hero />
                <Capabilities />
                <FirstAnalysis />
                <p className={styles.terms}>
                    GPL-3.0-or-later. Distribution Statement A: approved for public release,
                    distribution unlimited. <Link to="/public-release">Copyright and License</Link>
                </p>
            </main>
        </Layout>
    );
}
