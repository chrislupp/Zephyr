import type { ReactNode } from "react";
import {
    ellipticalChords,
    feasibility,
    initialChords,
    optimality,
    optimizedChords,
    problem,
    stations,
} from "@site/src/data/chord-optimization";
import styles from "./ChordFigures.module.css";

// Figures 2 and 3, drawn from the committed output of the OpenMDAO driver in
// `examples/chord_optimization.py`, so they are the solved result rather than a
// drawing of it. Figure 2 shows the starboard half only: the wing is mirrored
// about y = 0 and the other half carries the same numbers.

const HALF_SPAN = problem.span / 2;

const round = (value: number): number => Math.round(value * 100) / 100;

// ── Figure 2: planform ───────────────────────────────────────────────────────
// Spanwise across, chordwise down, one scale for both so a planform looks
// like a planform. The quarter-chord line is held straight and unswept by the
// example, so it is the horizontal datum of the drawing.

const PLAN = {
    left: 96,
    top: 30,
    width: 540,
    /** Drawing units per unit length, shared by both axes. */
    scale: 540 / HALF_SPAN,
};

const planX = (y: number): number => PLAN.left + y * PLAN.scale;
/** Chordwise station, measured aft from the quarter-chord line. */
const planY = (x: number): number => PLAN.top + (x + 0.25 * maxChord) * PLAN.scale;

const maxChord = Math.max(...optimizedChords, ...ellipticalChords, ...initialChords);

/** A closed half-planform outline: leading edge out, trailing edge back. */
const planform = (chords: number[]): string => {
    const leading = chords.map((c, j) => `${round(planX(stations[j]))},${round(planY(-0.25 * c))}`);
    const trailing = chords
        .map((c, j) => `${round(planX(stations[j]))},${round(planY(0.75 * c))}`)
        .reverse();
    return [...leading, ...trailing].join(" ") + " " + leading[0];
};

/** The ellipse drawn as a curve rather than through the design stations. */
const ellipse = (): string => {
    const steps = 160;
    const root = (4 * problem.area) / (Math.PI * problem.span);
    const chord = (y: number): number => root * Math.sqrt(Math.max(0, 1 - (y / HALF_SPAN) ** 2));
    const leading: string[] = [];
    const trailing: string[] = [];
    for (let i = 0; i <= steps; i += 1) {
        const y = (HALF_SPAN * i) / steps;
        leading.push(`${round(planX(y))},${round(planY(-0.25 * chord(y)))}`);
        trailing.push(`${round(planX(y))},${round(planY(0.75 * chord(y)))}`);
    }
    return [...leading, ...trailing.reverse()].join(" ");
};

const PLAN_BOTTOM = planY(0.75 * maxChord) + 6;
/** Chordwise ticks, positive aft of the quarter-chord line. */
const CHORD_TICKS = [-0.5, 0, 0.5, 1];
const CHORD_AXIS_X = PLAN.left - 12;
const PLAN_KEY = PLAN_BOTTOM + 66;
const PLAN_HEIGHT = PLAN_KEY + 46;

export function PlanformFigure(): ReactNode {
    return (
        <svg
            className={styles.figure}
            viewBox={`0 0 ${PLAN.left + PLAN.width + 24} ${PLAN_HEIGHT}`}
            role="img"
            aria-label="Plan view of the starboard half wing: the rectangular starting planform, the optimized planform, and the elliptical planform of the same span and area, which the optimized outline follows except at the tip."
        >
            <polyline className={styles.start} points={planform(initialChords)} />
            <polyline className={styles.reference} points={ellipse()} />
            <polyline className={styles.result} points={planform(optimizedChords)} />

            {optimizedChords.map((c, j) => (
                <circle
                    key={stations[j]}
                    className={styles.node}
                    cx={planX(stations[j])}
                    cy={planY(-0.25 * c)}
                    r={2.2}
                />
            ))}

            {/* span axis */}
            <line
                className={styles.axis}
                x1={PLAN.left}
                y1={PLAN_BOTTOM}
                x2={planX(HALF_SPAN)}
                y2={PLAN_BOTTOM}
            />
            {[0, 1, 2, 3, 4, 5].map((y) => (
                <g key={y}>
                    <line
                        className={styles.tick}
                        x1={planX(y)}
                        y1={PLAN_BOTTOM}
                        x2={planX(y)}
                        y2={PLAN_BOTTOM + 5}
                    />
                    <text
                        className={styles.label}
                        x={planX(y)}
                        y={PLAN_BOTTOM + 20}
                        textAnchor="middle"
                    >
                        {y}
                    </text>
                </g>
            ))}
            <text
                className={styles.label}
                x={PLAN.left + PLAN.width / 2}
                y={PLAN_BOTTOM + 44}
                textAnchor="middle"
            >
                Spanwise position (m)
            </text>

            {/* chordwise axis, positive aft, at the same scale as the span */}
            <line
                className={styles.axis}
                x1={CHORD_AXIS_X}
                y1={planY(CHORD_TICKS[0])}
                x2={CHORD_AXIS_X}
                y2={planY(CHORD_TICKS[CHORD_TICKS.length - 1])}
            />
            {CHORD_TICKS.map((x) => (
                <g key={x}>
                    <line
                        className={styles.tick}
                        x1={CHORD_AXIS_X - 5}
                        y1={planY(x)}
                        x2={CHORD_AXIS_X}
                        y2={planY(x)}
                    />
                    <text
                        className={styles.label}
                        x={CHORD_AXIS_X - 9}
                        y={planY(x) + 4}
                        textAnchor="end"
                    >
                        {x.toFixed(1).replace("-", "\u2212")}
                    </text>
                </g>
            ))}
            <text
                className={styles.label}
                x={18}
                y={(planY(CHORD_TICKS[0]) + planY(CHORD_TICKS[CHORD_TICKS.length - 1])) / 2}
                textAnchor="middle"
                transform={`rotate(-90 18 ${(planY(CHORD_TICKS[0]) + planY(CHORD_TICKS[CHORD_TICKS.length - 1])) / 2})`}
            >
                Chordwise position (m)
            </text>

            <Key y={PLAN_KEY} x={PLAN.left} />
        </svg>
    );
}

// ── Figure 3: IPOPT convergence ────────────────────────────────────────────
// IPOPT's own convergence measures, read from its output file: dual
// infeasibility (optimality) above, primal infeasibility (feasibility) below,
// each on a logarithmic axis, against iteration.

const PANEL = { left: 96, width: 500, height: 150 };
const TOP = 12;
const GAP = 64;
const LOWER = TOP + PANEL.height + GAP;
const LAST = feasibility.length - 1;

type Axis = { top: number; low: number; high: number; ticks: number[] };

/** Decade range of a panel, from 10^low to 10^high. */
const OPTIMALITY_AXIS: Axis = { top: TOP, low: -9, high: -1, ticks: [-1, -3, -5, -7, -9] };
const FEASIBILITY_AXIS: Axis = { top: LOWER, low: -16, high: -1, ticks: [-1, -4, -7, -10, -13, -16] };

const iterX = (index: number): number => PANEL.left + (index / LAST) * PANEL.width;

const logY = (axis: Axis, value: number): number => {
    const exponent = Math.min(axis.high, Math.max(axis.low, Math.log10(Math.max(value, 1e-300))));
    return axis.top + ((axis.high - exponent) / (axis.high - axis.low)) * PANEL.height;
};

const CONVERGENCE_HEIGHT = LOWER + PANEL.height + 74;

function LogPanel({ axis, values, label }: { axis: Axis; values: number[]; label: string }): ReactNode {
    const points = values.map((value, index) => [iterX(index), logY(axis, value)]);
    return (
        <g>
            {axis.ticks.map((exponent) => (
                <g key={exponent}>
                    <line
                        className={styles.tick}
                        x1={PANEL.left - 5}
                        y1={logY(axis, 10 ** exponent)}
                        x2={PANEL.left}
                        y2={logY(axis, 10 ** exponent)}
                    />
                    <text
                        className={styles.label}
                        x={PANEL.left - 10}
                        y={logY(axis, 10 ** exponent) + 4}
                        textAnchor="end"
                    >
                        10
                        <tspan dy={-6} fontSize="10">
                            {`\u2212${-exponent}`}
                        </tspan>
                    </text>
                </g>
            ))}
            <line
                className={styles.axis}
                x1={PANEL.left}
                y1={axis.top}
                x2={PANEL.left}
                y2={axis.top + PANEL.height}
            />
            <line
                className={styles.axis}
                x1={PANEL.left}
                y1={axis.top + PANEL.height}
                x2={PANEL.left + PANEL.width}
                y2={axis.top + PANEL.height}
            />
            <polyline
                className={styles.result}
                points={points.map(([x, y]) => `${round(x)},${round(y)}`).join(" ")}
            />
            {points.map(([x, y], index) => (
                <circle key={index} className={styles.node} cx={x} cy={y} r={3} />
            ))}
            <text
                className={styles.label}
                x={20}
                y={axis.top + PANEL.height / 2}
                textAnchor="middle"
                transform={`rotate(-90 20 ${axis.top + PANEL.height / 2})`}
            >
                {label}
            </text>
        </g>
    );
}

export function ConvergenceFigure(): ReactNode {
    const step = Math.max(5, 5 * Math.ceil(LAST / 30));
    const iterations = Array.from(
        { length: Math.floor(LAST / step) + 1 },
        (_, index) => index * step,
    );
    const bottom = LOWER + PANEL.height;
    return (
        <svg
            className={styles.figure}
            viewBox={`0 0 ${PANEL.left + PANEL.width + 24} ${CONVERGENCE_HEIGHT}`}
            role="img"
            aria-label="Two stacked logarithmic plots against IPOPT iteration: above, optimality; below, feasibility. Both fall below the tolerance by the final iteration."
        >
            <LogPanel axis={OPTIMALITY_AXIS} values={optimality} label="Optimality" />
            <LogPanel axis={FEASIBILITY_AXIS} values={feasibility} label="Feasibility" />
            {iterations.map((index) => (
                <g key={index}>
                    <line
                        className={styles.tick}
                        x1={iterX(index)}
                        y1={bottom}
                        x2={iterX(index)}
                        y2={bottom + 5}
                    />
                    <text
                        className={styles.label}
                        x={iterX(index)}
                        y={bottom + 20}
                        textAnchor="middle"
                    >
                        {index}
                    </text>
                </g>
            ))}
            <text
                className={styles.label}
                x={PANEL.left + PANEL.width / 2}
                y={bottom + 44}
                textAnchor="middle"
            >
                Iteration
            </text>
        </svg>
    );
}

// ── Key for the planform figure ─────────────────────────────────────────────

function Key({ x, y }: { x: number; y: number }): ReactNode {
    const entries: [string, string][] = [
        [styles.start, "rectangular start"],
        [styles.result, "optimized"],
        [styles.reference, "elliptical"],
    ];
    return (
        <g>
            {entries.map(([className, label], index) => {
                const row = x + index * 184;
                return (
                    <g key={label}>
                        <line className={className} x1={row} y1={y} x2={row + 30} y2={y} />
                        <text className={styles.label} x={row + 38} y={y + 4}>
                            {label}
                        </text>
                    </g>
                );
            })}
        </g>
    );
}
