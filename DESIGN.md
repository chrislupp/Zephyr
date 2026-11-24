---
name: Zephyr
description: The documentation as a standards manual. Every page is a numbered sheet on one grid.
colors:
  red: "#d8301f"
  red-text: "#c4281a"
  on-red: "#ffffff"
  select: "#f6cdc7"
  sheet: "#ffffff"
  rail: "#eef0f1"
  rail-press: "#e2e5e8"
  ink: "#14171a"
  ink-2: "#33383e"
  mute: "#5f666f"
  rule: "#c9cdd2"
  hair: "#e1e4e7"
  listing: "#f1f2f3"
  listing-highlight: "#e3e6e9"
  listing-string: "#1e6650"
  listing-name: "#22507f"
  listing-literal: "#7a4a00"
  listing-punctuation: "#4a5159"
  red-dark: "#ff5c47"
  red-text-dark: "#ff7361"
  on-red-dark: "#121518"
  select-dark: "#7a2a21"
  sheet-dark: "#121518"
  rail-dark: "#1a1e22"
  rail-press-dark: "#242a2f"
  ink-dark: "#e8eaec"
  ink-2-dark: "#c3c8ce"
  mute-dark: "#9aa1a9"
  rule-dark: "#3a4047"
  hair-dark: "#272c31"
  listing-dark: "#1a1e22"
  listing-highlight-dark: "#2a3036"
  listing-keyword-dark: "#ffffff"
  listing-string-dark: "#7fd1b0"
  listing-name-dark: "#8fbfea"
  listing-literal-dark: "#e9b872"
  listing-punctuation-dark: "#b4bac1"
typography:
  sheet-number:
    fontFamily: "Public Sans, Helvetica Neue, Helvetica, Arial, sans-serif"
    fontSize: "4.75rem"
    fontWeight: 300
    lineHeight: 1
    letterSpacing: "-0.03em"
    fontFeature: "tnum"
  display:
    fontFamily: "Public Sans, Helvetica Neue, Helvetica, Arial, sans-serif"
    fontSize: "3.75rem"
    fontWeight: 300
    lineHeight: 1.05
    letterSpacing: "-0.025em"
  display-api:
    fontFamily: "Public Sans, Helvetica Neue, Helvetica, Arial, sans-serif"
    fontSize: "3rem"
    fontWeight: 300
    lineHeight: 1.05
    letterSpacing: "-0.025em"
  headline:
    fontFamily: "Public Sans, Helvetica Neue, Helvetica, Arial, sans-serif"
    fontSize: "1.5rem"
    fontWeight: 600
    lineHeight: 1.2
    letterSpacing: "-0.01em"
  title:
    fontFamily: "Public Sans, Helvetica Neue, Helvetica, Arial, sans-serif"
    fontSize: "1.1875rem"
    fontWeight: 600
    lineHeight: 1.2
  subtitle:
    fontFamily: "Public Sans, Helvetica Neue, Helvetica, Arial, sans-serif"
    fontSize: "1rem"
    fontWeight: 600
    lineHeight: 1.2
  masthead:
    fontFamily: "Public Sans, Helvetica Neue, Helvetica, Arial, sans-serif"
    fontSize: "1.5rem"
    fontWeight: 700
    letterSpacing: "-0.02em"
  lede:
    fontFamily: "Public Sans, Helvetica Neue, Helvetica, Arial, sans-serif"
    fontSize: "1.0625rem"
    fontWeight: 400
    lineHeight: 1.6
  body:
    fontFamily: "Public Sans, Helvetica Neue, Helvetica, Arial, sans-serif"
    fontSize: "1rem"
    fontWeight: 400
    lineHeight: 1.6
    fontFeature: "tnum"
  rail:
    fontFamily: "Public Sans, Helvetica Neue, Helvetica, Arial, sans-serif"
    fontSize: "0.9375rem"
    fontWeight: 400
    lineHeight: 1.35
    fontFeature: "tnum"
  small:
    fontFamily: "Public Sans, Helvetica Neue, Helvetica, Arial, sans-serif"
    fontSize: "0.875rem"
    fontWeight: 400
    lineHeight: 1.5
  label:
    fontFamily: "Public Sans, Helvetica Neue, Helvetica, Arial, sans-serif"
    fontSize: "0.8125rem"
    fontWeight: 600
    lineHeight: 1.4
  code:
    fontFamily: "Source Code Pro, ui-monospace, SF Mono, Menlo, Consolas, monospace"
    fontSize: "0.875rem"
    fontWeight: 400
    lineHeight: 1.6
rounded:
  none: "0"
spacing:
  inset: "1.75rem"
  gutter: "4rem"
  gutter-closed: "1.5rem"
  number-col: "4.25rem"
  number-col-lettered: "4.875rem"
  measure: "38rem"
  wide: "46rem"
  rail-width: "20rem"
  contents-width: "15.5rem"
  masthead-height: "4.25rem"
components:
  masthead:
    backgroundColor: "{colors.sheet}"
    textColor: "{colors.ink}"
    typography: "{typography.masthead}"
    height: "{spacing.masthead-height}"
  brand-bar:
    backgroundColor: "{colors.red}"
    width: "0.75rem"
    height: "2.125rem"
  rail:
    backgroundColor: "{colors.rail}"
    width: "{spacing.rail-width}"
  rail-entry:
    textColor: "{colors.ink-2}"
    typography: "{typography.rail}"
    rounded: "{rounded.none}"
  rail-entry-hover:
    backgroundColor: "{colors.rail-press}"
    textColor: "{colors.ink}"
  rail-entry-current:
    backgroundColor: "{colors.sheet}"
    textColor: "{colors.ink}"
  place-mark:
    backgroundColor: "{colors.red}"
    textColor: "{colors.on-red}"
    width: "3.5rem"
  sheet-number:
    textColor: "{colors.red}"
    typography: "{typography.sheet-number}"
  sheet-title:
    textColor: "{colors.ink}"
    typography: "{typography.display}"
    width: "{spacing.wide}"
  sheet-title-api:
    textColor: "{colors.ink}"
    typography: "{typography.display-api}"
    width: "{spacing.wide}"
  section-number:
    textColor: "{colors.mute}"
    typography: "{typography.body}"
  listing:
    backgroundColor: "{colors.listing}"
    textColor: "{colors.ink}"
    typography: "{typography.code}"
    rounded: "{rounded.none}"
    padding: "1rem"
  listing-button:
    backgroundColor: "{colors.sheet}"
    textColor: "{colors.ink}"
    rounded: "{rounded.none}"
  mode-toggle:
    textColor: "{colors.ink}"
    rounded: "{rounded.none}"
  mode-toggle-hover:
    backgroundColor: "{colors.rail}"
  key-cap:
    backgroundColor: "{colors.listing}"
    textColor: "{colors.ink}"
    rounded: "{rounded.none}"
  note:
    textColor: "{colors.ink}"
    typography: "{typography.body}"
    rounded: "{rounded.none}"
    padding: "0.75rem 0"
  table:
    textColor: "{colors.ink}"
    typography: "{typography.rail}"
    padding: "0.5rem 1.5rem 0.5rem 0"
  pager-link:
    textColor: "{colors.ink}"
    rounded: "{rounded.none}"
    padding: "0.75rem 0 0"
  pager-number:
    textColor: "{colors.red-text}"
  contents-link:
    textColor: "{colors.ink-2}"
    typography: "{typography.small}"
  contents-mark:
    backgroundColor: "{colors.red}"
    size: "0.5rem"
  route-number:
    textColor: "{colors.red-text}"
  dropdown-menu:
    backgroundColor: "{colors.sheet}"
    textColor: "{colors.ink-2}"
    rounded: "{rounded.none}"
    padding: "0.25rem 0"
---

# Design System: Zephyr

## Overview

**Creative North Star: "The Standards Manual"**

The documentation is built as a 1970s federal graphics standards manual. Every page is a numbered sheet on one grid, and the sheet number is how a reader finds, cites, and returns to it. A bright white sheet sits beside a cool grey rail, the ink is near black, and one standards red is spent on numbers and place marks. One grotesque, Public Sans, sets everything except code and mathematics.

The frame commits hard and the reading column stays calm. That split is a decision, not an accident: the maintainer called the first, quieter sketches bland and then chose this direction, so the masthead rule, the rail, the red place mark, and the large sheet number are loud on purpose, while prose, listings, equations, and tables are set plainly at a fixed measure for an engineer reading in the middle of work. Light is the default sheet. Dark is the same sheet after hours, token for token, and is never designed separately.

This is a research library's reference. It does not read as a marketing site, as the stock Docusaurus template, or as the near-black, glowing-cyan theme it replaced.

**Key Characteristics:**
- Every sheet is numbered, and the rail, the head of the sheet, and the pager print the same number.
- One datum edge runs down the sheet; text registers to it and numbers hang in the gutter to its left.
- Structure comes from rules, not boxes. Corners are square and nothing casts a shadow.
- One red, only on sheet numbers, place marks, the brand bar, and the data ink of Plate 1. Hover and focus are said in ink.
- One grotesque in two registers: light at display size, regular to bold for text, with tabular numerals throughout.
- Motion answers state changes only; Plate 1's wake is the single authored movement.

## Colors

White sheet, grey rail, black ink, and one red; the only other hues are three quiet inks confined to code listings.

### Primary
- **Standards Red** (`red`, #d8301f): the fill and display red. It sets the sheet number above each title, fills the current sheet's number cell in the rail, the square beside the section in view, the rule under the current book in the masthead, the brand bar, and the bound vortices of Plate 1. It reaches 4.81:1 on the sheet, so it is used as text only at display size; small red text takes Standards Red Text.
- **Standards Red Text** (`red-text`, #c4281a): the same red darkened to stay legible as text, 5.73:1 on the sheet and 5.01:1 on the rail. It sets sheet numbers at text size (pager, route lists) and the numeric labels of the bound vortices on Plate 1.
- **On Red** (`on-red`, #ffffff): the number inside a filled place mark.
- **Selection Tint** (`select`, #f6cdc7): the red let down to a tint, for selected text only.

### Secondary
Listing inks. They appear inside code listings and nowhere else, and none of them is red.
- **String Green** (`listing-string`, #1e6650): strings, characters, attribute values.
- **Name Blue** (`listing-name`, #22507f): functions, types, macros, built-ins.
- **Literal Amber** (`listing-literal`, #7a4a00): numbers, booleans, constants, lifetimes.
- **Punctuation Grey** (`listing-punctuation`, #4a5159): punctuation, operators, namespaces.
- Keywords are set in Ink at weight 600, not in a colour. Comments are Mute and upright.

### Neutral
- **Sheet White** (`sheet`, #ffffff): the page, the masthead, the colophon, and the row of the current sheet in the rail.
- **Rail Grey** (`rail`, #eef0f1): the rail, the brand cell above it, and the hover ground of dropdown items, tab items, and the colour-mode toggle.
- **Rail Press** (`rail-press`, #e2e5e8): a rail row under the pointer.
- **Ink** (`ink`, #14171a): text, headings, links at rest and under the pointer, heavy rules, focus outlines, the text caret, and the page-load bar.
- **Second Ink** (`ink-2`, #33383e): rail entries, contents entries, captions, and quoted text.
- **Mute** (`mute`, #5f666f): hanging numbers, the running head, labels, equation numbers, and the tagline. 5.81:1 on the sheet, 5.08:1 on the rail.
- **Rule** (`rule`, #c9cdd2): light rules under group labels, the borders of the copy button and of a key cap, panel edges on a plate, and the scrollbar.
- **Hair** (`hair`, #e1e4e7): the datum line, table rows, and the closing rule of a block.
- **Listing Ground** (`listing`, #f1f2f3): the ground of code listings, inline code, and key caps. A highlighted line takes Listing Highlight (`listing-highlight`, #e3e6e9).

### Dark sheet
The dark variant redefines the same tokens under `[data-theme="dark"]`; the frontmatter carries each as its `-dark` sibling. The sheet becomes #121518, the rail and the listing ground #1a1e22, ink #e8eaec, and the red lifts to #ff5c47 (text #ff7361) so it holds 6.0:1 on the dark sheet. On Red inverts to the sheet colour, so the number in a place mark is dark on red. Listing keywords become pure white at weight 600; the three listing inks lift to #7fd1b0, #8fbfea, and #e9b872.

### Named Rules
**The Standards Red Rule.** Red appears only on sheet numbers (above the title, in the pager, in the route lists), place marks (the rail's number cell, the square beside the section in view, the rule under the current book), the brand bar, and the keyed bound vortices of Plate 1. It is never a syntax colour, never a note or warning colour, never a hover or focus colour, and never decoration. The theme's own primary variable (`--ifm-color-primary`) is set to ink, so any stock state the build has not restyled falls back to ink rather than red.

**The Never Colour Alone Rule.** A place mark always changes something besides hue: the current rail row also turns sheet white and bold, the current contents entry and the current book turn bold. Hover follows the same thought from the other side: it changes weight of line, not hue.

**The Same Sheet Rule.** Light is the default and dark is its variant. Every colour is a `--zm-*` token redefined once for dark; no component carries a colour of its own for either theme, except the listing inks, which are set as a pair of themes in the site config.

## Typography

**Display Font:** Public Sans (with Helvetica Neue, Helvetica, Arial, sans-serif)
**Body Font:** Public Sans (same stack)
**Label/Mono Font:** Source Code Pro (with ui-monospace, SF Mono, Menlo, Consolas, monospace)

**Character:** One government grotesque does the whole job, light and large for numbers and titles, plain and sturdy for text, so the page reads as an issued standard and not as a brand. Source Code Pro is its drafting-room companion and is used for code alone.

Both families are self-hosted variable fonts (Latin and Latin Extended; Public Sans upright and italic, Source Code Pro upright). No italic of the code face ships, so nothing in a listing is italic, comments included. The root size is 106.25%, so 1rem is 17px and every size below is in those rems. Tabular numerals (`"tnum"`) are on for the whole document. Code ligatures are off. Mathematics is set by KaTeX in its own faces at 1.1em.

### Hierarchy
- **Sheet number** (300, 4.75rem, line-height 1, -0.03em, Standards Red): the number above every sheet title. 3rem below 997px.
- **Display** (300, 3.75rem, line-height 1.05, -0.025em, balanced wrap): the sheet title and the name on the contents page. 2.5rem below 997px (2.75rem for the contents-page name). A title runs to the 46rem sheet width, not the prose measure, and never breaks inside a word.
- **Display, API sheets** (300, 3rem): API sheets are titled with identifiers, which must not break, so they are set a step smaller: 2.25rem below 997px, and 1.5rem (-0.015em) at 640px and below.
- **Headline** (600, 1.5rem, line-height 1.2, -0.01em): numbered section headings, and the route and Conventions headings on the contents page.
- **Title** (600, 1.1875rem, line-height 1.2): subsections. A fourth level sits at 600, 1rem.
- **Masthead** (700, 1.5rem, -0.02em): the name beside the brand bar.
- **Body** (400, 1rem, line-height 1.6): prose, held to a 38rem measure, about 66 characters. The contents-page statement is set one step up (1.0625rem). Strong text is 700.
- **Rail** (400, 0.9375rem, line-height 1.35): rail entries, masthead links (500), tables, and tab labels. Section heads and the current entry are 700.
- **Small** (400, 0.875rem): on-page contents, captions, the colophon, the tagline, dropdown items.
- **Label** (600, 0.8125rem, Mute, sentence case): group labels such as the book heads and "On this page". The running head and listing language labels use the same size at 400.
- **Code** (Source Code Pro 400, 0.875rem, line-height 1.6): listings. Inline code is 0.875em on the listing ground with no border.

### Named Rules
**The One Grotesque Rule.** Public Sans sets every word that is not code or mathematics, including labels on listings and plates. Source Code Pro never sets prose, navigation, or labels.

**The Light At Size Rule.** Weight 300 is reserved for the sheet number, the sheet title, and the contents-page name. Nothing else is lighter than 400. The only light type below 2.25rem is an API sheet's identifier title on a sheet 640px or narrower, which shrinks so that it never breaks.

**The Sentence Case Rule.** Labels are sentence case at normal tracking. The one exception in the build is a note's heading, which is 700, uppercase, tracked 0.04em.

## Layout

One grid, read from the left: rail, gutter, datum, column, contents.

The rail is 20rem wide on the rail ground and holds a number column (4.25rem, or 4.875rem in the lettered API book) that starts at a 1.75rem inset. Right of the rail a 4rem gutter leads to the datum, a single 1px hair line down the whole sheet. Content begins 1.75rem right of the datum, and that edge is shared by the running head, the sheet number, the title, headings, prose, listings, equations, tables, and the pager. Section numbers hang in the gutter left of the datum, right-aligned. The masthead (4.25rem, closed by a 2px ink rule) puts the name over the rail on the rail's ground and the tagline at the gutter; the colophon sets its one line on the same content edge.

On a sheet the column is at most 46rem wide. Prose, section headings, notes, and equations hold a 38rem measure inside it; the sheet title, listings, and tables may run the full 46rem but start on the same edge. Vertical rhythm is set by the headings (2.75rem above a section heading and 0.75rem below; 2rem and 0.5rem for a subsection), 1.125rem between paragraphs, and 1.5rem above and 1.75rem below a listing, table, or note.

The grid gives ground in steps:
- **1440px and up:** a 15.5rem on-page contents column sits right of the text. On the contents page the statement and Plate 1 share a row.
- **Below 1440px:** the on-page contents fold above the text as a ruled, collapsible block, because the 38rem measure would not survive beside them.
- **Below 1280px:** the gutter closes to 1.5rem and section numbers move in front of their headings, with the heading's second line returning under its own first letter. The tagline leaves the masthead below 1200px.
- **Below 997px:** the rail leaves the sheet (a drawer on doc sheets; below the text on the contents page and the distribution statement), side padding is 1rem, and the sheet number and title step down.
- **640px and below:** Plate 1 switches to the starboard half of the wing at nearly twice the scale, the pager's two neighbours stack so neither title is cut, and API sheet titles take their smallest size. At 480px and below display equations lose their indent.

### Named Rules
**The Datum Rule.** Every block on a sheet starts on the one edge right of the datum line. Nothing is centred, and the only indents from it are list items, a display equation (2rem), and quoted text (1.75rem).

**The 66-Character Rule.** Prose holds 38rem at every width. A wider element is a sheet title, a listing, a table, or a plate, never a paragraph.

**The No-Zero Rule.** Sections and sheets count from 1: sections are 1, 2, 3 and their sheets 1.1, 1.2. The API reference book is lettered, A1, A2.1. The contents page carries no sheet number.

## Elevation & Depth

The system is flat. Every shadow token is set to none, and the build adds no shadow anywhere. Depth is carried by ground alone: the grey rail beside the white sheet, the listing ground under code, and the current rail row turning sheet white, as if the sheet showed through the rail at the reader's place. The few things that float over the sheet (the version dropdown, a listing's copy button, the back-to-top control) separate with a 1px border, not a shadow. Stock parts are brought onto the sheet the same way: a key cap is flat, and the page-load bar is a plain ink line with no glow.

### Named Rules
**The Flat Sheet Rule.** No shadows and no glows, at rest or in any state. Separation is a change of ground or a rule.

## Shapes

Every corner is square (radius 0, set globally and again on each component that had its own). The form language is the ruled line in four weights:

- **Heavy rule** (2px, Ink): opens a block. It closes the masthead, tops every listing, table, note, plate, route, and pager link, and opens the colophon.
- **Ink hairline** (1px, Ink): the rule under a table head, and the full border of a floating control (dropdown menu, back-to-top).
- **Light rule** (1px, Rule): under group labels and listing titles, and around the copy button and a key cap.
- **Hair** (1px, Hair): the datum, table rows, and the closing rule of a note or a plate.

Marks are rectangles: the brand bar (0.75rem by 2.125rem), the filled number cell in the rail, the 0.5rem square beside the section in view, and the 4px rule under the current book. The favicon is the same vocabulary: a white square, a red bar, and an ink Z. The colour-mode toggle is square. Circles appear only inside a plate, as data. The focus outline is 2px Ink at a 2px offset.

### Named Rules
**The Ruled Block Rule.** A block opens with a 2px ink rule and, if it closes, closes with a hair. Content is never enclosed on four sides; a full border belongs only to a small control floating over the sheet.

## Components

Hover is said in ink. A link's underline thickens from 1px to 2px, a rail row takes a darker ground, a pager link's rule thickens; nothing turns red. State changes that ease do so at 160ms on `cubic-bezier(0.16, 1, 0.3, 1)` and touch only colour, background, opacity, and the rail caret's rotation. Reduced motion collapses every transition and animation. There is no call-to-action button and no form field in the system; the entries below are everything the build has.

### Buttons
The only buttons are small controls that float over the sheet.
- **Shape:** square (radius 0), with a 1px border.
- **Listing copy button:** Sheet White ground, Ink glyph, 1px Rule border. Hidden until the listing is hovered or the button takes focus.
- **Back to top:** Sheet White ground, 1px Ink border, Ink glyph.
- **Colour-mode toggle:** the stock control made square, Ink glyph, no border.
- **Hover / Focus:** the copy button's border turns Ink; back-to-top and the colour-mode toggle take the Rail Grey ground. Focus is the 2px Ink outline.

### Cards / Containers
There are no cards. Three ruled blocks do the containing.
- **Listing:** Listing Ground, a 2px ink rule on top, no border and no radius, 1rem padding. A label in the text face (0.8125rem, Mute, "Rust", "Shell", "TOML") names the language above the code; a listing title, when given, is 600 with a light rule beneath.
- **Note:** no ground at all. A 2px ink rule above, a hair below, 0.75rem vertical padding, held to the 38rem measure. The heading is 0.8125rem, 700, uppercase. Every kind of note, including warnings and the unreleased-version banner, is set in Ink; kinds are not colour-coded.
- **Table:** a 2px ink rule above the head and a 1px ink rule below it, hair between rows, no vertical rules, no stripes, no head ground. Cells are left-aligned and top-aligned at 0.9375rem with 0.5rem by 1.5rem padding; head cells are 600 and do not wrap. A wide table scrolls sideways inside itself.
- **Definition rows** (contents page): the same ruling applied to a term and its definition, a 6.5rem term column in 600.
- **Key cap:** Listing Ground, a 1px light rule, square, no shadow, set in the code face in Ink.

### Navigation
- **Masthead:** Sheet White, 4.25rem, closed by a 2px ink rule. The brand is the red bar and the name; from 997px its cell is rail-wide and takes the rail ground on any page with a rail. Links sit at the right (0.9375rem, 500, Ink) and take a 2px ink underline at a 0.35em offset under the pointer. The current book is 700 with a 4px red rule laid over the masthead's own rule. The version dropdown is Mute at 400 and opens a square menu with a 1px ink border.
- **Rail:** every entry is a number cell and a label. Section heads are 700 Ink and carry the section number alone; sheets are 400 Second Ink and carry section and sheet. Hover fills the row with Rail Press. Section heads fold their sheets with a drawn chevron; every section starts open, so the rail reads as a contents page.
- **Place mark (rail):** the current sheet's number cell fills Standards Red with the number in On Red, and the row turns Sheet White and 700.
- **Running head:** one line above the sheet number, the path at the datum in Mute with "/" separators and the current step in Ink, the version label at the column's far edge. A step under the pointer turns Ink and underlines. The home step is omitted because the masthead already says it.
- **On-page contents:** a labelled list at 0.875rem, numbered to match the section numbers on the sheet; sub-entries are Mute. An entry under the pointer turns Ink and underlines. The entry in view is 700 Ink with a 0.5rem red square in the margin beside it.
- **Pager:** two links, each under a 2px ink rule, with a Mute sublabel and a 600 title preceded by the neighbouring sheet's number in Standards Red Text. On hover the rule thickens to 4px and the text holds its place; the title stays Ink. Titles wrap between words only. On a sheet 640px or narrower the two links stack.
- **Links in text:** Ink, underlined 1px at a 0.2em offset; the underline thickens to 2px under the pointer. The colophon carries no links: it is the copyright line under the masthead's twin, a 2px ink rule, and nothing else.

### Sheet Number
The signature. A light, red number (4.75rem, 300) stands above the sheet title on the same edge. It is computed from the sheet's position in its sidebar, the same list the rail counts with CSS counters, so on a doc sheet the two agree by construction, and the pager prints it again for the neighbouring sheets. A screen reader hears "Sheet" before it. The distribution statement sits outside the versioned docs but is still a sheet: it takes the rail, its number, and the place mark, and its rail starts at the same height as on a doc sheet. Its number, and the numbers in the contents page's route lists, are looked up from the site's own contents list and never typed in.

### Hanging Numbers
Section headings on a sheet are counted and the number is printed in Mute at 400: in the gutter left of the datum from 1280px, in front of the heading below that. Display equations are counted the same way and carry "(n)" in Mute at the right edge of the measure, vertically centred on the equation, which is left-aligned at a 2rem indent and scrolls sideways if it is too long.

### Route List
The contents page's two ways in. Each route opens with a 2px ink rule, a headline, and one line of Second Ink, then lists its sheets as hair-ruled rows: the sheet number in Standards Red Text in a 3rem column, the title at 600, underlined in ink under the pointer. Two columns from 997px. These are rules, not cards.

### Plate
A figure is a plate. It opens with a 2px ink rule and closes with a hair, and its caption (0.875rem, Second Ink) begins with "Plate n." in bold Ink. Plate 1 is drawn from the solver's own mesh for the Quick Start wing, not illustrated: geometry in Ink (1px, 1.5px for the outline), panel edges in Rule, bound vortices in Standards Red (1.5px, 3px in the key), control points as Ink dots, labels in the text face in Mute. Strokes do not scale with the drawing; the trailing legs are the exception, because their dashes are measured in drawing units so that each leg is drawn whole at any size. Its trailing legs draw aft once on load, root first (900ms each on the house ease, 28ms apart); with reduced motion they are simply drawn. Only this one plate exists; a second should follow it.

## Do's and Don'ts

### Do:
- **Do** number every new sheet by its place in the rail, counting from 1, and print that same number in the rail, above the title, and in the pager.
- **Do** start every block on the content edge right of the datum, hold prose to 38rem, and let only listings, tables, and plates run wider.
- **Do** open a block with a 2px ink rule and close it with a 1px hair; keep every corner square.
- **Do** pair every red place mark with a second cue: a change of weight or of ground.
- **Do** say hover in ink: thicken the underline or the rule, or darken the ground.
- **Do** take every colour from a `--zm-*` token so the dark sheet follows without a second design.
- **Do** set keywords in heavier ink and keep listings to green strings, blue names, and amber literals.
- **Do** label things in the text face, in sentence case, in Mute: listing languages, group heads, plate annotations.
- **Do** keep state changes to 160ms on `cubic-bezier(0.16, 1, 0.3, 1)` and honour reduced motion.

### Don't:
- **Don't** put red on anything but sheet numbers, place marks, the brand bar, and a plate's data ink. It stays out of syntax colours, notes, warnings, and hover and focus states.
- **Don't** start any numbering at 0, and don't give the contents page a sheet number.
- **Don't** round a corner, cast a shadow, or add a glow.
- **Don't** enclose content in a four-sided box or a card. A full border belongs only to a small control floating over the sheet.
- **Don't** colour-code notes by kind. A warning is ruled and set in Ink like any other note.
- **Don't** add a second text face, and don't set prose, navigation, or labels in Source Code Pro.
- **Don't** use weight 300 for anything but the sheet number and the title.
- **Don't** break a title inside a word; set an identifier title smaller instead.
- **Don't** italicise anything in a listing; the code face ships upright only.
- **Don't** make dark the default or give it colours of its own outside the tokens.
- **Don't** animate anything that is not a change of state. Plate 1's wake is the one authored movement.
- **Don't** let a page read as marketing: no slogan hero, no feature cards, and nothing of the stock template or the earlier near-black, glowing-cyan theme.
