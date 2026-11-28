import type { PrismTheme } from "prism-react-renderer";
import type { Config } from "@docusaurus/types";
import type * as Preset from "@docusaurus/preset-classic";
import remarkMath from "remark-math";
import rehypeKatex from "rehype-katex";

// Listing colours. Dark is the default theme; light is its variant.
const listingLight: PrismTheme = {
    plain: { color: "#0D1520", backgroundColor: "#FFFFFF" },
    styles: [
        { types: ["comment", "prolog", "doctype", "cdata"], style: { color: "#5F6B78" } },
        { types: ["keyword", "atrule", "important", "deleted"], style: { color: "#6B2FA3" } },
        { types: ["string", "char", "attr-value", "inserted"], style: { color: "#0B6B4F" } },
        { types: ["function", "class-name", "builtin", "macro", "tag"], style: { color: "#1F55B5" } },
        { types: ["number", "boolean", "constant", "symbol", "lifetime"], style: { color: "#9A4A00" } },
        { types: ["punctuation", "operator", "namespace"], style: { color: "#46525F", opacity: 1 } },
    ],
};

const listingDark: PrismTheme = {
    plain: { color: "#E6EDF3", backgroundColor: "#10161D" },
    styles: [
        { types: ["comment", "prolog", "doctype", "cdata"], style: { color: "#7D8C9B" } },
        { types: ["keyword", "atrule", "important", "deleted"], style: { color: "#C4A7F5" } },
        { types: ["string", "char", "attr-value", "inserted"], style: { color: "#6FDCB6" } },
        { types: ["function", "class-name", "builtin", "macro", "tag"], style: { color: "#8FC3FF" } },
        { types: ["number", "boolean", "constant", "symbol", "lifetime"], style: { color: "#F2C46D" } },
        { types: ["punctuation", "operator", "namespace"], style: { color: "#A9B6C3", opacity: 1 } },
    ],
};

const config: Config = {
    title: "Zephyr",
    tagline: "Vortex lattice method library",
    favicon: "img/favicon.svg",

    url: "https://chrislupp.github.io",
    baseUrl: "/Zephyr/",

    organizationName: "chrislupp",
    projectName: "Zephyr",

    onBrokenLinks: "throw",

    i18n: {
        defaultLocale: "en",
        locales: ["en"],
    },

    presets: [
        [
            "classic",
            {
                docs: {
                    sidebarPath: "./sidebars.ts",
                    // The live docs are the default version, served at
                    // /docs/. Released snapshots keep their own numbered
                    // paths (/docs/0.5.0/...) and stay reachable from the
                    // version dropdown. Without this, /docs/ would resolve to
                    // the newest frozen snapshot, which cannot carry a page
                    // added since that release.
                    lastVersion: "current",
                    // Every page starts listed: the rail is a contents page.
                    sidebarCollapsed: false,
                    remarkPlugins: [remarkMath],
                    rehypePlugins: [rehypeKatex],
                    editUrl:
                        "https://github.com/chrislupp/Zephyr/tree/develop/docs/",
                    versions: {
                        current: {
                            // The live docs track the development branch; released
                            // snapshots carry their own numbers in the dropdown.
                            label: "Development",
                        },
                    },
                },
                blog: false,
                theme: {
                    // KaTeX ships from node_modules so the build carries the
                    // stylesheet that matches the renderer; no CDN, no hash to
                    // keep in step.
                    customCss: [
                        require.resolve("katex/dist/katex.min.css"),
                        "./src/css/custom.css",
                    ],
                },
            } satisfies Preset.Options,
        ],
    ],

    themeConfig: {
        colorMode: {
            // Dark is the default, as in a post-processor's render window;
            // the switch in the navbar offers the light variant.
            defaultMode: "dark",
            respectPrefersColorScheme: false,
        },
        navbar: {
            title: "Zephyr",
            items: [
                {
                    type: "docSidebar",
                    sidebarId: "docsSidebar",
                    position: "right",
                    label: "Docs",
                },
                {
                    type: "docSidebar",
                    sidebarId: "apiSidebar",
                    position: "right",
                    label: "API",
                },
                {
                    to: "/docs/examples/basic-analysis",
                    label: "Examples",
                    position: "right",
                    activeBaseRegex: "/docs/examples/",
                },
                {
                    type: "docsVersionDropdown",
                    position: "right",
                },
                {
                    href: "https://github.com/chrislupp/Zephyr",
                    label: "GitHub",
                    position: "right",
                },
            ],
        },
        footer: {
            // The footer is the copyright line and nothing else. Copyright,
            // licence terms, and the distribution statement live on the
            // Copyright and License page, which the rail lists on every doc
            // page and the start page links to.
            copyright: `Copyright © 2024-${new Date().getFullYear()} Christopher A. Lupp.`,
        },
        prism: {
            theme: listingLight,
            darkTheme: listingDark,
            additionalLanguages: ["rust", "bash", "toml"],
        },
    } satisfies Preset.ThemeConfig,
};

export default config;
