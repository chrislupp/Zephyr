import type { CSSProperties, ReactNode } from "react";
import styles from "./LatticePlate.module.css";

// Plate 1 is drawn from the lattice Zephyr builds for the Quick Start wing:
// 4 chordwise by 16 spanwise panels, default cosine spacing, span 10, chord 1.
// The stations below are that mesh's own values (AVL CSPACER chordwise, sine
// of equal angle steps spanwise), so the drawing is the wing as solved, not an
// illustration of one.
const PANELS_SPAN = 16;
const HALF_SPAN = 5;
const CHORD_EDGES = [0, 0.178606, 0.5, 0.821394, 1];
const BOUND_VORTEX = [0.030154, 0.25, 0.586824, 0.883022];
const CONTROL_POINT = [0.116978, 0.413176, 0.75, 0.969846];

const spanStation = (k: number): number =>
    HALF_SPAN * Math.sin(((k - PANELS_SPAN / 2) * Math.PI) / PANELS_SPAN);

const ALL_EDGES = Array.from({ length: PANELS_SPAN + 1 }, (_, k) => ({
    k,
    y: spanStation(k),
}));
const ALL_CENTRES = Array.from({ length: PANELS_SPAN }, (_, k) => spanStation(k + 0.5));

const LEFT = 24;
const LEADING_EDGE = 44;
const LEG_SOLID = 86;
const LEG_FADE = 30;

const round = (value: number): number => Math.round(value * 100) / 100;

type Frame = {
    /** Width of the plan and of the chordwise detail, in SVG units. */
    width: number;
    /** Spanwise station at the plan's left edge: -5 for the whole wing, 0 for the starboard half. */
    from: number;
    viewBox: string;
    spanLabel: string;
    titleLines: string[];
    label: string;
};

// The whole wing for a wide sheet; the starboard half, at nearly twice the
// scale, for a narrow one. Both keep the plan and the detail on one width.
const WIDE: Frame = {
    width: 560,
    from: -HALF_SPAN,
    viewBox: "0 0 664 460",
    spanLabel: "span 10",
    titleLines: ["One strip, leading edge to trailing edge, in fractions of chord"],
    label: "Plan view of the Quick Start wing divided into 4 chordwise by 16 spanwise panels, with a bound vortex and a control point in each panel and trailing legs running aft; below it, the chordwise stations of one strip.",
};

const NARROW: Frame = {
    width: 300,
    from: 0,
    viewBox: "0 0 400 484",
    spanLabel: "half span 5",
    titleLines: ["One strip, leading edge to trailing edge,", "in fractions of chord"],
    label: "Plan view of the starboard half of the Quick Start wing, 4 chordwise by 8 spanwise panels, with a bound vortex and a control point in each panel and trailing legs running aft; below it, the chordwise stations of one strip.",
};

function Plate({ frame, className }: { frame: Frame; className: string }): ReactNode {
    const { width, from } = frame;
    const right = LEFT + width;
    const scale = width / (HALF_SPAN - from);
    const trailingEdge = LEADING_EDGE + scale;
    const legEnd = trailingEdge + LEG_SOLID + LEG_FADE;
    const axesX = LEFT + 12;
    const axesY = legEnd + 34;
    const titleY = axesY + 86;
    const detailY = titleY + 54 + (frame.titleLines.length - 1) * 18;
    const legendY = detailY + 56;
    const rootEdge = PANELS_SPAN / 2;

    const planX = (y: number): number => LEFT + (y - from) * scale;
    const planY = (xOverC: number): number => LEADING_EDGE + xOverC * scale;
    const detailX = (xOverC: number): number => LEFT + xOverC * width;

    const edges = ALL_EDGES.filter(({ y }) => y >= from - 1e-9);
    const centres = ALL_CENTRES.filter((y) => y > from);
    // Edges strictly inside the outline; the outline draws the outer two.
    const inner = edges.filter(({ y }) => planX(y) > LEFT + 0.5 && planX(y) < right - 0.5);

    return (
        <svg
            className={`${styles.plate} ${className}`}
            viewBox={frame.viewBox}
            role="img"
            aria-label={frame.label}
        >
            {/* Span and chord dimensions. */}
            <g className={styles.ink}>
                <line x1={LEFT} y1="26" x2={right} y2="26" />
                <line x1={LEFT} y1="20" x2={LEFT} y2="32" />
                <line x1={right} y1="20" x2={right} y2="32" />
                <line x1={right + 14} y1={LEADING_EDGE} x2={right + 14} y2={trailingEdge} />
                <line x1={right + 8} y1={LEADING_EDGE} x2={right + 20} y2={LEADING_EDGE} />
                <line x1={right + 8} y1={trailingEdge} x2={right + 20} y2={trailingEdge} />
            </g>
            <text x={LEFT + width / 2} y="15" textAnchor="middle" className={styles.label}>
                {frame.spanLabel}
            </text>
            <text x={right + 24} y={LEADING_EDGE + scale / 2 + 4} className={styles.label}>
                chord 1
            </text>

            {/* Trailing legs, one from every strip edge, running aft. */}
            <g className={styles.legs}>
                {edges.map(({ y, k }) => (
                    <line
                        key={`leg-${k}`}
                        className={styles.leg}
                        style={{ "--zm-order": Math.abs(k - rootEdge) } as CSSProperties}
                        x1={round(planX(y))}
                        y1={trailingEdge}
                        x2={round(planX(y))}
                        y2={trailingEdge + LEG_SOLID}
                    />
                ))}
                {edges.map(({ y, k }) => (
                    <line
                        key={`fade-${k}`}
                        className={styles.legFade}
                        x1={round(planX(y))}
                        y1={trailingEdge + LEG_SOLID}
                        x2={round(planX(y))}
                        y2={legEnd}
                    />
                ))}
            </g>
            <text x={right} y={legEnd + 20} textAnchor="end" className={styles.label}>
                trailing legs continue aft
            </text>

            {/* Body axes: x aft, y starboard, z up (out of the sheet). */}
            <g className={styles.ink}>
                <line x1={axesX} y1={axesY} x2={axesX} y2={axesY + 44} />
                <path
                    d={`M${axesX - 4} ${axesY + 38} L${axesX} ${axesY + 46} L${axesX + 4} ${axesY + 38}`}
                    fill="none"
                />
                <line x1={axesX} y1={axesY} x2={axesX + 44} y2={axesY} />
                <path
                    d={`M${axesX + 38} ${axesY - 4} L${axesX + 46} ${axesY} L${axesX + 38} ${axesY + 4}`}
                    fill="none"
                />
                <circle cx={axesX} cy={axesY} r="5" className={styles.sheet} />
            </g>
            <circle cx={axesX} cy={axesY} r="1.5" className={styles.fill} />
            <text x={axesX + 54} y={axesY + 4} className={styles.label}>
                y starboard
            </text>
            <text x={axesX + 10} y={axesY + 46} className={styles.label}>
                x aft
            </text>
            <text x={axesX - 8} y={axesY - 12} className={styles.label}>
                z up
            </text>

            {/* Panel edges. */}
            <g className={styles.hair}>
                {inner.map(({ y, k }) => (
                    <line
                        key={`span-${k}`}
                        x1={round(planX(y))}
                        y1={LEADING_EDGE}
                        x2={round(planX(y))}
                        y2={trailingEdge}
                    />
                ))}
                {CHORD_EDGES.slice(1, -1).map((x) => (
                    <line
                        key={`chord-${x}`}
                        x1={LEFT}
                        y1={round(planY(x))}
                        x2={right}
                        y2={round(planY(x))}
                    />
                ))}
            </g>
            <rect
                className={styles.outline}
                x={LEFT}
                y={LEADING_EDGE}
                width={width}
                height={scale}
            />

            {/* Bound vortices and control points. */}
            <g className={styles.red}>
                {BOUND_VORTEX.map((x) => (
                    <line
                        key={`bound-${x}`}
                        x1={LEFT}
                        y1={round(planY(x))}
                        x2={right}
                        y2={round(planY(x))}
                    />
                ))}
            </g>
            <g className={styles.fill}>
                {centres.flatMap((y) =>
                    CONTROL_POINT.map((x) => (
                        <circle
                            key={`cp-${y}-${x}`}
                            cx={round(planX(y))}
                            cy={round(planY(x))}
                            r="1.6"
                        />
                    )),
                )}
            </g>

            {/* One strip, leading edge to trailing edge. */}
            <text x={LEFT} y={titleY} className={styles.title}>
                {frame.titleLines.map((line, i) => (
                    <tspan key={line} x={LEFT} dy={i === 0 ? 0 : 18}>
                        {line}
                    </tspan>
                ))}
            </text>
            <line className={styles.outline} x1={LEFT} y1={detailY} x2={right} y2={detailY} />
            <g className={styles.ink}>
                {CHORD_EDGES.map((x) => (
                    <line
                        key={`tick-${x}`}
                        x1={round(detailX(x))}
                        y1={detailY - 7}
                        x2={round(detailX(x))}
                        y2={detailY + 7}
                    />
                ))}
            </g>
            <g className={styles.redBar}>
                {BOUND_VORTEX.map((x) => (
                    <line
                        key={`bar-${x}`}
                        x1={round(detailX(x))}
                        y1={detailY - 14}
                        x2={round(detailX(x))}
                        y2={detailY + 14}
                    />
                ))}
            </g>
            <g className={styles.fill}>
                {CONTROL_POINT.map((x) => (
                    <circle key={`dot-${x}`} cx={round(detailX(x))} cy={detailY} r="4" />
                ))}
            </g>
            {BOUND_VORTEX.map((x, i) => (
                <text
                    key={`bv-${x}`}
                    x={round(detailX(x))}
                    y={detailY - 24}
                    textAnchor={i === 0 ? "start" : "middle"}
                    className={styles.redLabel}
                >
                    {x.toFixed(3)}
                </text>
            ))}
            {CONTROL_POINT.map((x, i) => (
                <text
                    key={`cpl-${x}`}
                    x={round(detailX(x))}
                    y={detailY + 34}
                    textAnchor={i === CONTROL_POINT.length - 1 ? "end" : "middle"}
                    className={styles.label}
                >
                    {x.toFixed(3)}
                </text>
            ))}

            {/* Key. */}
            <g className={styles.redBar}>
                <line x1={LEFT + 2} y1={legendY - 11} x2={LEFT + 2} y2={legendY + 3} />
            </g>
            <text x={LEFT + 12} y={legendY} className={styles.label}>
                bound vortex
            </text>
            <circle cx={LEFT + 112} cy={legendY - 4} r="4" className={styles.fill} />
            <text x={LEFT + 122} y={legendY} className={styles.label}>
                control point
            </text>
            <g className={styles.ink}>
                <line x1={LEFT + 222} y1={legendY - 11} x2={LEFT + 222} y2={legendY + 3} />
            </g>
            <text x={LEFT + 232} y={legendY} className={styles.label}>
                panel edge
            </text>
        </svg>
    );
}

/** Both compositions ship; the stylesheet shows the one that fits the sheet. */
export default function LatticePlate(): ReactNode {
    return (
        <>
            <Plate frame={WIDE} className={styles.wide} />
            <Plate frame={NARROW} className={styles.narrow} />
        </>
    );
}
