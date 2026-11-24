import type { ReactNode } from "react";
import clsx from "clsx";
import Link from "@docusaurus/Link";
import type { Props } from "@theme/PaginatorNavLink";
import { useSheetNumberFor } from "@site/src/components/SheetNumber";

// The stock pager link, with the neighbouring sheet's number printed beside
// its title so "what comes next" is a sheet a reader can cite.
export default function PaginatorNavLink(props: Props): ReactNode {
    const { permalink, title, subLabel, isNext } = props;
    const number = useSheetNumberFor(permalink);
    return (
        <Link
            className={clsx(
                "pagination-nav__link",
                isNext ? "pagination-nav__link--next" : "pagination-nav__link--prev",
            )}
            to={permalink}
        >
            {subLabel && <div className="pagination-nav__sublabel">{subLabel}</div>}
            <div className="pagination-nav__label">
                {number && <span className="zm-pager-number">{number}</span>}
                {title}
            </div>
        </Link>
    );
}
