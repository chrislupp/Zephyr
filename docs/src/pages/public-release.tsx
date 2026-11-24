import type { ReactNode } from "react";
import Head from "@docusaurus/Head";
import Layout from "@theme/Layout";
import MDXContent from "@theme/MDXContent";
import ContentsRail from "@site/src/components/ContentsRail";
import { docs, sheetNumber } from "@site/src/data/contents";
import Statement, * as statement from "./_public-release.md";

// Copyright, license, and the distribution statement sit outside the versioned
// docs so they apply to every version, but this is still a sheet of the manual:
// it takes the rail, its number, and the place mark. The footer carries the
// copyright line alone and links here for the rest. The text lives in
// `_public-release.md`.
// The MDX loader exports the file's front matter; the ambient `*.md` type
// only declares the component.
const { frontMatter } = statement as unknown as {
    frontMatter: { title: string; description: string };
};

const PATH = "/public-release";
const number = sheetNumber(docs, PATH);

export default function PublicRelease(): ReactNode {
    return (
        <Layout title={frontMatter.title} description={frontMatter.description}>
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
