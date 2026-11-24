import type { ReactNode } from "react";
import Content from "@theme-original/DocItem/Content";
import type ContentType from "@theme/DocItem/Content";
import type { WrapperProps } from "@docusaurus/types";
import SheetNumber from "@site/src/components/SheetNumber";

type Props = WrapperProps<typeof ContentType>;

// Prints the sheet number above the page title; everything else is the stock
// component.
export default function ContentWrapper(props: Props): ReactNode {
    return (
        <>
            <SheetNumber />
            <Content {...props} />
        </>
    );
}
