import type { ReactNode } from "react";
import { useDocsSidebar } from "@docusaurus/plugin-content-docs/client";
import { useLocation } from "@docusaurus/router";

const trimSlash = (path: string): string => path.replace(/\/+$/, "") || "/";

/**
 * The number of the current page in its sidebar: "2" for a top-level entry,
 * "2.1" for the first page of the second section. The rail prints the same
 * numbers with CSS counters over the same list, so the two always agree. API
 * reference sheets take an "A" prefix, matching `--zm-sheet-prefix` in
 * custom.css.
 */
export function useSheetNumberFor(path: string): string | null {
    const sidebar = useDocsSidebar();
    const here = trimSlash(path);
    if (!sidebar) {
        return null;
    }
    const prefix = sidebar.name === "apiSidebar" ? "A" : "";
    let section = 0;
    for (const item of sidebar.items) {
        section += 1;
        if (item.type === "html") {
            continue;
        }
        if (item.href && trimSlash(item.href) === here) {
            return `${prefix}${section}`;
        }
        if (item.type !== "category") {
            continue;
        }
        let sheet = 0;
        for (const child of item.items) {
            sheet += 1;
            if (
                child.type !== "html" &&
                child.href &&
                trimSlash(child.href) === here
            ) {
                return `${prefix}${section}.${sheet}`;
            }
        }
    }
    return null;
}

export default function SheetNumber(): ReactNode {
    const number = useSheetNumberFor(useLocation().pathname);
    if (!number) {
        return null;
    }
    return (
        <div className="zm-sheetno">
            <span className="zm-sr-only">Sheet </span>
            {number}
        </div>
    );
}
