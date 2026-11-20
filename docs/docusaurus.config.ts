import { themes as prismThemes } from "prism-react-renderer";
import type { Config } from "@docusaurus/types";
import type * as Preset from "@docusaurus/preset-classic";
import remarkMath from "remark-math";
import rehypeKatex from "rehype-katex";

const config: Config = {
    title: "VLM",
    tagline:
        "Vortex Lattice Method Library for Aerodynamic Analysis",
    favicon: "img/favicon.ico",

    url: "https://chrislupp.github.io",
    baseUrl: "/Zephyr/",

    organizationName: "chrislupp",
    projectName: "Zephyr",

    onBrokenLinks: "throw",

    i18n: {
        defaultLocale: "en",
        locales: ["en"],
    },

    stylesheets: [
        {
            href: "https://cdn.jsdelivr.net/npm/katex@0.16.11/dist/katex.min.css",
            type: "text/css",
            integrity:
                "sha384-nB0miv6/jRmo5RLHO8BIp/8hwC1slNFDuv3VUgI4A7O/lEKmpRgiuNOQI2bRpFB",
            crossorigin: "anonymous",
        },
    ],

    presets: [
        [
            "classic",
            {
                docs: {
                    sidebarPath: "./sidebars.ts",
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
                    customCss: "./src/css/custom.css",
                },
            } satisfies Preset.Options,
        ],
    ],

    themeConfig: {
        colorMode: {
            defaultMode: "dark",
            respectPrefersColorScheme: true,
        },
        navbar: {
            title: "VLM",
            items: [
                {
                    type: "docSidebar",
                    sidebarId: "docsSidebar",
                    position: "left",
                    label: "Docs",
                },
                {
                    type: "docSidebar",
                    sidebarId: "apiSidebar",
                    position: "left",
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
            style: "dark",
            links: [
                {
                    title: "Documentation",
                    items: [
                        {
                            label: "Getting Started",
                            to: "/docs/getting-started/installation",
                        },
                        {
                            label: "API Reference",
                            to: "/docs/api",
                        },
                    ],
                },
                {
                    title: "More",
                    items: [
                        {
                            label: "GitHub",
                            href: "https://github.com/chrislupp/Zephyr",
                        },
                        {
                            label: "Changelog",
                            href: "https://github.com/chrislupp/Zephyr/blob/main/CHANGELOG.md",
                        },
                        {
                            label: "Distribution Statement",
                            to: "/public-release",
                        },
                    ],
                },
            ],
            copyright: `Copyright \u00A9 2024-${new Date().getFullYear()} Christopher A. Lupp. Approved for public release; distribution unlimited (AFRL-2026-2198). Built with Docusaurus.`,
        },
        prism: {
            theme: prismThemes.github,
            darkTheme: prismThemes.dracula,
            additionalLanguages: ["rust", "bash"],
        },
    } satisfies Preset.ThemeConfig,
};

export default config;
