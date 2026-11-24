import type { CSSProperties, ReactNode } from "react";
import Link from "@docusaurus/Link";
import { api, docs, type Section } from "@site/src/data/contents";
import styles from "./ContentsRail.module.css";

// The rail for pages that live outside the docs plugin (the contents page and
// the distribution statement). It uses the docs sidebar's own class names, so
// the theme's rail rules and its counters number it exactly as a book's rail.

const sameSheet = (a: string, b?: string): boolean => b !== undefined && a === b;

function Book({
    sections,
    prefix,
    current,
}: {
    sections: Section[];
    prefix: string;
    current?: string;
}): ReactNode {
    return (
        <ul
            className="theme-doc-sidebar-menu menu__list"
            style={
                {
                    "--zm-sheet-prefix": `"${prefix}"`,
                    "--zm-number-col": prefix ? "4.875rem" : "4.25rem",
                } as CSSProperties
            }
        >
            {sections.map((section) => (
                <li
                    key={section.label}
                    className="theme-doc-sidebar-item-link menu__list-item"
                >
                    {section.to ? (
                        <Link
                            className={`menu__link${sameSheet(section.to, current) ? " menu__link--active" : ""}`}
                            aria-current={sameSheet(section.to, current) ? "page" : undefined}
                            to={section.to}
                        >
                            <span>{section.label}</span>
                        </Link>
                    ) : (
                        <div className="menu__list-item-collapsible">
                            <span className={`menu__link ${styles.sectionHead}`}>
                                <span>{section.label}</span>
                            </span>
                        </div>
                    )}
                    {section.sheets.length > 0 && (
                        <ul className="menu__list">
                            {section.sheets.map((sheet) => (
                                <li
                                    key={sheet.to}
                                    className="theme-doc-sidebar-item-link menu__list-item"
                                >
                                    <Link className="menu__link" to={sheet.to}>
                                        <span>{sheet.label}</span>
                                    </Link>
                                </li>
                            ))}
                        </ul>
                    )}
                </li>
            ))}
        </ul>
    );
}

export default function ContentsRail({
    withApi = false,
    current,
}: {
    /** Also list the API reference book, under a heading for each book. */
    withApi?: boolean;
    /** The `to` of the sheet being shown, so its number cell is filled. */
    current?: string;
}): ReactNode {
    return (
        <aside className={`theme-doc-sidebar-container ${styles.rail}`}>
            <nav className="menu" aria-label="Contents">
                {withApi && <h2 className={styles.book}>Documentation</h2>}
                <Book sections={docs} prefix="" current={current} />
                {withApi && (
                    <>
                        <h2 className={styles.book}>API reference</h2>
                        <Book sections={api} prefix="A" current={current} />
                    </>
                )}
            </nav>
        </aside>
    );
}
