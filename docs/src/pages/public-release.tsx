import type { ReactNode } from "react";
import Layout from "@theme/Layout";
import MDXContent from "@theme/MDXContent";
import ContentsRail from "@site/src/components/ContentsRail";
import Statement from "./_public-release.md";

// Copyright, license, and the distribution statement sit outside the versioned
// docs so they apply to every version, but this is still a page of the
// documentation: it takes the same sidebar. The footer carries the copyright
// line alone; the sidebar lists this page for the rest. The text lives in
// `_public-release.md`. That file is an MDX partial, which Docusaurus forbids
// from carrying front matter (a hard error when CI is set), so the page
// metadata lives here.
const TITLE = "Copyright and License";
const DESCRIPTION =
    "Copyright, license, and the public release and distribution statement for Zephyr.";

const PATH = "/public-release";

export default function PublicRelease(): ReactNode {
    return (
        <Layout title={TITLE} description={DESCRIPTION}>
            <div className="z-page">
                <ContentsRail current={PATH} />
                <main className="z-page__main">
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
