import type { PrismTheme } from "prism-react-renderer";
import type { Config } from "@docusaurus/types";
import type * as Preset from "@docusaurus/preset-classic";
import remarkMath from "remark-math";
import rehypeKatex from "rehype-katex";

// Listing colours for the two sheets. Red is reserved for sheet numbers and
// place marks, so keywords are set in heavier ink; green is for strings, blue
// for names, amber for literals.
const listingLight: PrismTheme = {
    plain: { color: "#14171A", backgroundColor: "#F1F2F3" },
    styles: [
        {
            types: ["comment", "prolog", "doctype", "cdata"],
            style: { color: "#5F666F" },
        },
        {
            types: ["keyword", "atrule", "important", "deleted"],
            style: { color: "#14171A", fontWeight: "600" },
        },
        {
            types: ["string", "char", "attr-value", "inserted"],
            style: { color: "#1E6650" },
        },
        {
            types: ["function", "class-name", "builtin", "macro", "tag"],
            style: { color: "#22507F" },
        },
        {
            types: ["number", "boolean", "constant", "symbol", "lifetime"],
            style: { color: "#7A4A00" },
        },
        {
            types: ["punctuation", "operator", "namespace"],
            style: { color: "#4A5159", opacity: 1 },
        },
    ],
};

const listingDark: PrismTheme = {
    plain: { color: "#E8EAEC", backgroundColor: "#1A1E22" },
    styles: [
        {
            types: ["comment", "prolog", "doctype", "cdata"],
            style: { color: "#9AA1A9" },
        },
        {
            types: ["keyword", "atrule", "important", "deleted"],
            style: { color: "#FFFFFF", fontWeight: "600" },
        },
        {
            types: ["string", "char", "attr-value", "inserted"],
            style: { color: "#7FD1B0" },
        },
        {
            types: ["function", "class-name", "builtin", "macro", "tag"],
            style: { color: "#8FBFEA" },
        },
        {
            types: ["number", "boolean", "constant", "symbol", "lifetime"],
            style: { color: "#E9B872" },
        },
        {
            types: ["punctuation", "operator", "namespace"],
            style: { color: "#B4BAC1", opacity: 1 },
        },
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
                    // The live docs are the default sheet set, served at
                    // /docs/. Released snapshots keep their own numbered
                    // paths (/docs/0.5.0/...) and stay reachable from the
                    // version dropdown. Without this, /docs/ would resolve to
                    // the newest frozen snapshot, which cannot carry a sheet
                    // added since that release.
                    lastVersion: "current",
                    // Every sheet starts listed: the rail is a contents page.
                    sidebarCollapsed: false,
                    remarkPlugins: [remarkMath],
                    rehypePlugins: [rehypeKatex],
                    editUrl:
                        "https://github.com/chrislupp/Zephyr/tree/develop/docs/",
                    versions: {
                        current: {
                            label: "Next",
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
            defaultMode: "light",
            respectPrefersColorScheme: true,
        },
        navbar: {
            title: "Zephyr",
            items: [
                {
                    type: "html",
                    position: "left",
                    value: '<span class="zm-tagline">Vortex lattice method library</span>',
                },
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
            // The colophon is the copyright line and nothing else. Copyright,
            // licence terms, and the distribution statement live on the
            // Copyright and License sheet, which the rail lists on every page.
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
