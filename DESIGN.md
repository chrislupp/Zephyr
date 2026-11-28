---
name: Zephyr
description: The documentation as a post-processor. A graphite result window, near-white ink, and the solution's viridis colour map as the only saturated colour.
colors:
  ground: "#0a0e13"
  ground-glow: "#16202b"
  surface: "#10161d"
  raised: "#161e27"
  ink: "#e6edf3"
  ink-2: "#a9b6c3"
  ink-3: "#8796a5"
  line: "#1f2a35"
  line-strong: "#2c3a48"
  accent: "#4fd1a5"
  accent-ink: "#6fdcb6"
  accent-soft: "rgba(79, 209, 165, 0.1)"
  select: "rgba(79, 209, 165, 0.26)"
  install-fill: "rgba(18, 26, 35, 0.8)"
  white: "#ffffff"
  ground-light: "#f6f8fa"
  surface-light: "#ffffff"
  raised-light: "#eef2f5"
  ink-light: "#0d1520"
  ink-2-light: "#46525f"
  ink-3-light: "#5f6b78"
  line-light: "#dde3e9"
  line-strong-light: "#c4ccd5"
  accent-light: "#0f7a69"
  accent-ink-light: "#0b5e51"
  accent-soft-light: "rgba(15, 122, 105, 0.1)"
  select-light: "rgba(33, 144, 141, 0.22)"
  viridis-00: "#440154"
  viridis-13: "#472c7a"
  viridis-25: "#3b518b"
  viridis-38: "#2c718e"
  viridis-50: "#21908d"
  viridis-63: "#27ad81"
  viridis-75: "#5cc863"
  viridis-88: "#aadc32"
  viridis-100: "#fde725"
  caution: "#b7791f"
  listing-comment: "#7d8c9b"
  listing-keyword: "#c4a7f5"
  listing-string: "#6fdcb6"
  listing-name: "#8fc3ff"
  listing-literal: "#f2c46d"
  listing-punctuation: "#a9b6c3"
  listing-comment-light: "#5f6b78"
  listing-keyword-light: "#6b2fa3"
  listing-string-light: "#0b6b4f"
  listing-name-light: "#1f55b5"
  listing-literal-light: "#9a4a00"
  listing-punctuation-light: "#46525f"
typography:
  display:
    fontFamily: "Host Grotesk, ui-sans-serif, system-ui, -apple-system, Segoe UI, sans-serif"
    fontSize: "clamp(2.5rem, 4.6vw, 4.25rem)"
    fontWeight: 620
    lineHeight: 1.02
    letterSpacing: "-0.035em"
  headline:
    fontFamily: "Host Grotesk, ui-sans-serif, system-ui, -apple-system, Segoe UI, sans-serif"
    fontSize: "2.75rem"
    fontWeight: 650
    lineHeight: 1.08
    letterSpacing: "-0.03em"
  section-title:
    fontFamily: "Host Grotesk, ui-sans-serif, system-ui, -apple-system, Segoe UI, sans-serif"
    fontSize: "clamp(1.75rem, 3vw, 2.5rem)"
    fontWeight: 620
    lineHeight: 1.1
    letterSpacing: "-0.03em"
  title:
    fontFamily: "Host Grotesk, ui-sans-serif, system-ui, -apple-system, Segoe UI, sans-serif"
    fontSize: "1.625rem"
    fontWeight: 650
    lineHeight: 1.2
    letterSpacing: "-0.02em"
  title-sm:
    fontFamily: "Host Grotesk, ui-sans-serif, system-ui, -apple-system, Segoe UI, sans-serif"
    fontSize: "1.25rem"
    fontWeight: 650
    lineHeight: 1.2
    letterSpacing: "-0.01em"
  lede:
    fontFamily: "Host Grotesk, ui-sans-serif, system-ui, -apple-system, Segoe UI, sans-serif"
    fontSize: "1.1875rem"
    fontWeight: 400
    lineHeight: 1.5
  body:
    fontFamily: "Host Grotesk, ui-sans-serif, system-ui, -apple-system, Segoe UI, sans-serif"
    fontSize: "1.0625rem"
    fontWeight: 400
    lineHeight: 1.65
  button:
    fontFamily: "Host Grotesk, ui-sans-serif, system-ui, -apple-system, Segoe UI, sans-serif"
    fontSize: "1rem"
    fontWeight: 600
  label:
    fontFamily: "Host Grotesk, ui-sans-serif, system-ui, -apple-system, Segoe UI, sans-serif"
    fontSize: "0.8125rem"
    fontWeight: 600
  caption:
    fontFamily: "Host Grotesk, ui-sans-serif, system-ui, -apple-system, Segoe UI, sans-serif"
    fontSize: "0.8125rem"
    fontWeight: 400
    lineHeight: 1.5
  readout:
    fontFamily: "Red Hat Mono, ui-monospace, SF Mono, Menlo, Consolas, monospace"
    fontSize: "1.25rem"
    fontWeight: 500
  code:
    fontFamily: "Red Hat Mono, ui-monospace, SF Mono, Menlo, Consolas, monospace"
    fontSize: "0.875rem"
    fontWeight: 400
    lineHeight: 1.65
rounded:
  xs: "4px"
  sm: "5px"
  md: "6px"
  base: "8px"
  control: "10px"
  frame: "12px"
  pill: "999px"
spacing:
  stack-xs: "0.375rem"
  stack-sm: "0.75rem"
  stack-md: "1.25rem"
  stack-lg: "1.5rem"
  stack-xl: "2.5rem"
  gap-xl: "3rem"
  gutter: "3.5rem"
  gutter-mobile: "1.125rem"
  section: "6rem"
  section-mobile: "4rem"
  measure: "42rem"
  rail: "18rem"
components:
  button-primary:
    backgroundColor: "{colors.ink}"
    textColor: "{colors.ground}"
    typography: "{typography.button}"
    rounded: "{rounded.control}"
    padding: "0 1.25rem"
    height: "3rem"
  button-primary-hover:
    backgroundColor: "{colors.white}"
    textColor: "{colors.ground}"
  button-install:
    backgroundColor: "{colors.install-fill}"
    textColor: "{colors.ink}"
    typography: "{typography.code}"
    rounded: "{rounded.control}"
    padding: "0 1.25rem"
    height: "3rem"
  button-secondary:
    textColor: "{colors.ink}"
    typography: "{typography.button}"
    rounded: "{rounded.control}"
    padding: "0 1.25rem"
    height: "3rem"
  menu-link-active:
    backgroundColor: "{colors.accent-soft}"
    textColor: "{colors.accent-ink}"
    rounded: "{rounded.md}"
    padding: "0.375rem 0.75rem"
  code-block:
    backgroundColor: "{colors.surface}"
    textColor: "{colors.ink}"
    typography: "{typography.code}"
    rounded: "{rounded.base}"
    padding: "1.125rem"
  inline-code:
    backgroundColor: "{colors.raised}"
    rounded: "{rounded.sm}"
    padding: "0.1em 0.35em"
  table-head:
    backgroundColor: "{colors.raised}"
    textColor: "{colors.ink-2}"
    typography: "{typography.label}"
    padding: "0.625rem 0.875rem"
  alert:
    backgroundColor: "{colors.surface}"
    textColor: "{colors.ink}"
    rounded: "{rounded.base}"
  pager-link:
    backgroundColor: "{colors.surface}"
    rounded: "{rounded.base}"
  figure-frame:
    backgroundColor: "{colors.surface}"
    rounded: "{rounded.base}"
    padding: "1rem 1rem 0.5rem"
  result-visual:
    backgroundColor: "{colors.ground}"
    rounded: "{rounded.frame}"
    padding: "1.25rem"
  version-badge:
    textColor: "{colors.ink-3}"
    rounded: "{rounded.pill}"
  colour-bar:
    rounded: "{rounded.xs}"
    height: "8px"
    width: "19rem"
  brand-mark:
    rounded: "{rounded.sm}"
    size: "1.25rem"
---

# Design System: Zephyr

## Overview

**Creative North Star: "The Post-processor"**

The site is a post-processor's result window. The ground is graphite, the ink is near-white, and the only saturated colour is the viridis colour map that Zephyr's own solutions are coloured with: the data's colour bar, from which the single accent, a viridis teal-green, is drawn for links and the reader's current place. Dark is the default (`colorMode.defaultMode: "dark"`, `respectPrefersColorScheme: false`); light is the variant, the same structure on a cool white ground, offered by the navbar switch.

The start page opens on one result: the uCRM transport wing, panelled and solved by Zephyr, filling a dark viewport with the offer laid over it, then three capabilities each shown by real output, then one listing. Doc pages keep the standard Docusaurus structure (rail, page, on-page contents) in the same palette. Density is a reference tool's: compact chrome, quiet 1px cool borders, prose held to a reading measure, and figures, tables and code allowed full width.

Depth comes from tone, not shadow. Surfaces step from ground to surface to raised; nothing glows. Motion is a 150ms colour or border change, plus one resolve on the start page where the solution fades in over its own panel mesh.

**Key Characteristics:**
- Graphite ground, near-white ink, viridis as the one saturated family.
- One accent, a viridis teal-green, for links, active rail entries, focus rings and the caret.
- Host Grotesk for text and headings; Red Hat Mono for code and numeric readouts.
- Small radii (4 to 12px), 1px cool borders, flat surfaces.
- Every image is the library's own output, rendered from a solved case.
- The start page's result window is dark in both themes.

## Colors

A graphite instrument palette: cool neutrals in two themes, one viridis-derived accent, and the full viridis map reserved for data.

### Primary
- **Viridis Teal** (`accent`): links (through `accent-ink`), the active rail and on-page contents entry, focus outlines, the caret, the hover border of pager links, and the solved planform on the start page. In light it becomes **Deep Viridis Teal** (`accent-light`, with `accent-ink-light` for link text so it holds contrast on white).
- **Teal Wash** (`accent-soft`): the fill behind the active rail entry and the active dropdown entry. **Teal Selection** (`select`) is the text-selection fill.

### Tertiary
- **The Viridis Map** (`viridis-00` to `viridis-100`): the solution's colour map, nine stops from deep purple to yellow. It appears as data in the renders, as the colour bar under the readouts, as the brand chip before the site name, and as the page-load progress bar. `viridis-100` also strokes the elliptic target planform on the start page.
- **Caution Amber** (`caution`): the border of warning and danger notes, and nothing else.

### Neutral
- **Graphite Ground** (`ground`) and **Graphite Glow** (`ground-glow`): the page ground, and the centre of the radial falloff in the result window and the capability visuals.
- **Panel Surface** (`surface`): code blocks, notes, pager links, figure frames, the mobile contents box.
- **Raised Panel** (`raised`): inline code, table heads, keyboard keys, menu hover.
- **Near-white Ink** (`ink`), **Secondary Ink** (`ink-2`), **Muted Ink** (`ink-3`): headings and body; rail and nav links, captions, blockquotes; breadcrumbs, footnotes, contents links, readout labels.
- **Cool Rule** (`line`) and **Strong Rule** (`line-strong`): every 1px border; the strong rule outlines buttons, keyboard keys and the GitHub link, and colours the scrollbar.
- Each neutral has a `-light` counterpart for the light theme: a cool white ground, white surfaces, blue-black ink.
- **Listing colours** (`listing-*`, `listing-*-light`): the Prism themes in `docusaurus.config.ts`, muted violet keywords, teal strings, blue names, gold literals, on the surface colour.

### Named Rules
**The Colour Bar Rule.** Outside code listings and the caution border, the only saturated colour on any page is viridis: the map itself on data, and the one teal accent drawn from it. A new hue needs a data reason.

**The Result Window Rule.** The start page's hero and the three capability visuals are dark in both themes. Their renders carry pale wake lines and viridis polygons on a transparent ground that only read on graphite, so the hero pins its own palette (`--h-ground`, `--h-ink`, and so on, equal to the dark tokens) and the navbar over it is forced to dark ink. Everything outside those frames follows the theme.

## Typography

**Display Font:** Host Grotesk (with ui-sans-serif, system-ui, -apple-system, Segoe UI, sans-serif)
**Body Font:** Host Grotesk
**Label/Mono Font:** Red Hat Mono (with ui-monospace, SF Mono, Menlo, Consolas, monospace)

**Character:** A contemporary engineering grotesk, set tight at display sizes and open in body, paired with a mono that reads as an instrument readout. Both are self-hosted variable faces (weights 300 to 800) in `docs/src/fonts/`.

### Hierarchy
- **Display** (`display`): the start page's headline only, held to 15ch and balanced over two lines.
- **Headline** (`headline`): doc page titles. API pages drop to 2.25rem and break anywhere so long item paths fit.
- **Section title** (`section-title`): the start page's section headings, held to 22ch.
- **Title** (`title`): doc page h2, with 3rem above. **Title small** (`title-sm`): doc h3 and capability names, with 2.25rem above.
- **Lede** (`lede`): the start page's one sentence under the headline, 33rem wide; 1.0625rem on phones.
- **Body** (`body`): doc prose at 1.0625rem, line height 1.65, held to the 42rem measure.
- **Label** (`label`): table heads, rail category heads, code block titles, admonition heads (with 0.02em tracking).
- **Caption** (`caption`): render captions, readout labels, breadcrumbs, footer, contents links. Figure captions are 0.875rem.
- **Readout** (`readout`): the start page's numeric results (C_L, C_Di, the validation reference).
- **Code** (`code`): listings at 0.875rem; inline code at 0.875em, ligatures off.

### Named Rules
**The Readout Rule.** Red Hat Mono is for code and for numbers a solver produced. It is never a label face or a decorative voice.

**The Tight Display Rule.** Tracking tightens as size grows: -0.035em at display, -0.03em at headline, -0.02em at h2, -0.01em at h3, none in body.

## Layout

Doc pages use the Docusaurus three-column frame: an 18rem rail (`rail`) with a 1px right rule, the page, and the on-page contents on the right. Prose, headings, blockquotes, display equations, notes and figure captions are capped at the 42rem measure (`measure`); tables, code blocks and figure frames may run wider and scroll horizontally inside themselves. The navbar is 4rem high and translucent (88% ground with a blurred, saturated backdrop) over a 1px bottom rule. The Copyright and License page reuses the same rail through a custom contents component, stacked under the page below 997px.

The start page is a full-height result window (`min-height: max(40rem, min(100svh, 58rem))`) laid out on a two-row grid: copy at the top left, readouts and colour bar along the bottom, the render absolutely placed across the right and centre behind them. Below it, sections sit in an 82rem column with 3.5rem gutters (`gutter`) and 6rem between sections (`section`): a three-column capability row with 1.5rem gaps, then a 5:7 split of text and listing.

Spacing steps that recur: 0.375rem, 0.75rem, 1.25rem, 1.5rem, 2.5rem and 3rem. At 1100px the render drops lower, capabilities go to two columns and the listing stacks; at 760px the hero becomes a single column with the render in flow (130% wide, bleeding both edges), gutters shrink to 1.125rem and capabilities stack.

## Elevation & Depth

The system is flat. Depth is tonal: ground, then surface, then raised, each separated by a 1px cool rule. Code blocks, notes, keyboard keys, the back-to-top button and the progress peg all have their shadows removed, and Infima's low shadow is `none`. The only shadows left are Infima's medium and top-level ones on floating menus. The result window gets its depth from a soft radial falloff (`ground-glow` at 62% 58% to `ground` at 70%), not from light effects.

### Shadow Vocabulary
- **Floating menu, light** (`box-shadow: 0 8px 24px rgba(13, 21, 32, 0.12)`): dropdowns and other Infima overlays.
- **Floating menu, dark** (`box-shadow: 0 12px 32px rgba(0, 0, 0, 0.45)`): the same in dark.

### Named Rules
**The No-Glow Rule.** Nothing emits light. No coloured shadows, no halos, no neon edges; a surface separates from its ground by tone and a 1px rule.

## Shapes

Small, consistent corners. The base radius is 8px (`base`), used by code blocks, tables, notes, pager links, figure frames, dropdowns, the GitHub link and the colour-mode toggle. Tighter steps go on smaller things: 6px on rail and dropdown links and breadcrumbs, 5px on inline code, keyboard keys and the brand chip, 4px on the focus ring and the colour bar. Start-page buttons take 10px (`control`) and the capability visuals 12px (`frame`). The version badge is the one pill. Borders are always 1px; blockquotes take a 2px strong-rule left edge. Tables are wrapped in a rounded, bordered block with row rules only, no column rules and no stripes.

## Components

### Buttons
Start-page actions, quiet and solid.
- **Shape:** gently rounded (10px), 3rem tall, 1.25rem side padding, 600 weight.
- **Primary:** near-white fill with graphite text inside the result window ("Get started"); hover goes to pure white. Below the window, the themed primary (`primaryAlt`) fills with the theme's ink and text in the theme's ground, hovering to 88% opacity.
- **Install:** the `cargo add zephyr-vlm` command as an outlined button in mono over a translucent graphite fill, with a ruled-off "Copy" state that reads "Copied" for 1.6s; hover lifts the border to muted ink.
- **Secondary:** transparent with a 1px strong-rule border; hover lifts the border to muted ink.
- **Motion:** background, border and colour transition over 150ms on `cubic-bezier(0.2, 0.8, 0.2, 1)`.

### Navigation
- **Navbar:** brand at left (a 1.25rem viridis chip, then "Zephyr" at 1.1875rem/650), links right at 0.9375rem/500 in secondary ink, hovering and active to ink. GitHub sits in an outlined 8px button with its external-link icon hidden. On the start page the navbar is transparent, unblurred and borderless over the result window, its colours pinned dark.
- **Rail:** 0.9375rem links in secondary ink with 6px corners; the active page is accent-ink on the teal wash at 600 weight; level-one categories are 0.8125rem/600 ink headings. All categories start expanded.
- **On-page contents:** 0.8125rem muted links that turn accent-ink when active or hovered.
- **Pager:** surface boxes with a 1px rule whose border turns accent on hover.

### Code
- **Blocks:** surface background, 1px rule, no shadow, 8px corners, 1.125rem padding, Red Hat Mono 0.875rem/1.65. Titles are sans 0.8125rem/600 in secondary ink over a rule. Highlighted lines take a faint teal (`rgba(79, 209, 165, 0.12)` dark, `rgba(33, 144, 141, 0.14)` light).
- **Inline:** raised fill, 5px corners, no border, ligatures off.

### Tables
A bordered, rounded block that fits its content and scrolls when it must. Heads are raised, 0.8125rem/600 secondary ink, no wrap; cells are 0.9375rem/1.5 with bottom rules only. A table whose head row is empty is treated as a list and its head is hidden.

### Notes
Surface background, 1px strong-rule border, no shadow, held to the measure. Warning and danger notes swap the border for caution amber. Links inside are accent-ink and underlined.

### Result Window (signature)
The start page's hero. Two stacked renders of the uCRM wing (680 panels, α = 2.5°): the panel mesh and, over it, the solution coloured by panel loading, which fades in once over 1100ms after 300ms (shown at once under reduced motion). The render is masked so its trailing legs dissolve aft (up the screen) and it fades in from the left, so nothing crosses the copy. A muted caption at top right names the case; it sits outside the masked figure. The version dropdown is hidden on the start page and shown on doc pages, where the live docs are labelled "Development". Along the bottom, readouts in mono with muted labels, and the colour bar: an 8px viridis ramp, 19rem wide, with 0 / 0.5 / ≥ 1.0 ticks.

### Capability Visual (signature)
A 4:3 frame, 12px corners, 1px rule, graphite radial fill in both themes, holding real output: the wing-and-tail render, the chord-optimization planforms (start dashed in muted ink, elliptic target dotted in `viridis-100`, IPOPT result solid teal over a teal wash), or the OptVL fixture list with teal dots. A caption naming the source sits under it, then the name, one sentence, and plain links.

### Figures on doc pages
A surface frame with a 1px rule and 8px corners holding an SVG plot drawn from the library's output. Axes in ink, ticks in strong rule, the starting shape in strong rule, the solved result in ink at the heaviest weight, the analytical reference dashed in accent. Labels are 13.5px with tabular figures. Captions are 0.875rem secondary ink under the frame, with bold runs in ink.

## Do's and Don'ts

### Do:
- **Do** take every colour from the `--z-*` tokens so both themes follow; hard-code hexes only inside the dark-in-both-themes result window and capability visuals.
- **Do** keep viridis for data, the colour bar, the brand chip and the progress bar, and the teal accent for links, the current place and focus.
- **Do** render every image from a solved Zephyr case with `docs/scripts/render_loading.py`, keeping its provenance comment, and name the case beside the image.
- **Do** separate surfaces with a 1px cool rule and a tonal step.
- **Do** hold prose, headings, notes and captions to the 42rem measure and let tables, code and figures run wider.
- **Do** keep hover and focus to 150ms colour and border changes, with a 2px accent outline offset 2px for focus.
- **Do** honour reduced motion: the resolve and every transition collapse to nothing.

### Don't:
- **Don't** add glow, coloured shadows, or shadows on code blocks, notes or keys.
- **Don't** introduce a second accent hue or use viridis as decoration away from data.
- **Don't** place an illustration, stock image or mock screenshot; images are the library's output.
- **Don't** let the result window or capability visuals follow the light theme.
- **Don't** set labels, headings or prose in Red Hat Mono.
- **Don't** pin a crate version in prose, install lines or examples.
