import type { ReactNode } from "react";
import Head from "@docusaurus/Head";
import Layout from "@theme/Layout";
import MDXContent from "@theme/MDXContent";
import ContentsRail from "@site/src/components/ContentsRail";
import { docs, sheetNumber } from "@site/src/data/contents";
import Statement from "./_public-release.md";

// Copyright, license, and the distribution statement sit outside the versioned
// docs so they apply to every version, but this is still a sheet of the manual:
// it takes the rail, its number, and the place mark. The footer carries the
// copyright line alone and links here for the rest. The text lives in
// `_public-release.md`. That file is an MDX partial, which Docusaurus forbids
// from carrying front matter (a hard error when CI is set), so the page
// metadata lives here.
const TITLE = "Copyright and License";
const DESCRIPTION =
    "Copyright, license, and the public release and distribution statement for Zephyr.";

const PATH = "/public-release";
const number = sheetNumber(docs, PATH);

export default function PublicRelease(): ReactNode {
    return (
        <Layout title={TITLE} description={DESCRIPTION}>
            <Head>
                <html className="zm-has-rail" />
            </Head>
            <div className="zm-page">
                <ContentsRail current={PATH} />
                <main className="zm-page__main">
                    {number && (
                        <div className="zm-sheetno">
                            <span className="zm-sr-only">Sheet </span>
                            {number}
                        </div>
                    )}
                    <div className="theme-doc-markdown markdown">
                        <MDXContent>
                            <Statement />
                        </MDXContent>
                    </div>
                </main>
            </div>
        </Layout>
    );
}
