---
version: 1
slug: "docs"
primary_target: "docs"
related_targets: []
---

# Surface brief: documentation site (`docs/`)

Scope: the whole Docusaurus site, every doc page plus the docs index at the site root.
Visitor mode: Read.

## Reader and job

Aerospace engineers and researchers, mid-work, with an aircraft model or an
optimizer in mind. Two jobs: run a first analysis, and wire Zephyr into an MDO
loop. Content: prose, Rust and shell listings, KaTeX mathematics, tables.

Constraints: page text and structure stay as written; released snapshots under
`docs/versioned_docs/` are untouched; the public-release clearance line stays in
the footer; no crate version in prose. Must not read as a marketing site, as
stock Docusaurus, or as the earlier near-black glowing-cyan theme. The
maintainer called the first sketches bland before choosing this direction, so
the frame commits hard; the reading column stays calm.

## Direction contract

THESIS: The documentation is a standards manual. Every page is a numbered sheet
on one grid, and the sheet number is how a reader finds, cites, and returns to
it. It refuses the dark template with a slogan hero and feature cards.

OWN-WORLD: Bright white sheet, cool grey rail, near-black ink, one standards red
spent only on sheet numbers and place marks (and the keyed bound vortices of
Plate 1); hover and focus are said in ink. One grotesque (Public Sans): light
for titles, regular and semibold for text; Source Code Pro for code alone.
Square corners, rules instead of boxes, tabular numerals, numbers hanging in a
gutter beside a single datum edge that every title, heading, listing, equation,
and table registers to. Display equations carry numbers at the column edge.

STORY: The reader learns what Zephyr is and which conventions it uses, takes one
of two routes (first analysis, MDO), and always knows which sheet they are on
and which comes next.

FIRST VIEWPORT: Doc page: masthead with a red bar and the name over the rail,
navigation at right, a 2px ink rule beneath. Rail lists every sheet with its
number; the current sheet's number cell is filled red. Column: running head,
the sheet number in red light grotesque, the title in light grotesque beneath
it, prose at a fixed 66-character measure; numbered contents at right where
the measure survives beside them. Contents page: no sheet number (the
maintainer ruled that nothing starts at 0); the name, the one-paragraph
statement, `cargo add zephyr-vlm`, Plate 1 (the Quick Start wing drawn as a
dimensioned lattice), the two routes, and the conventions; both books are
listed in its rail. Below 1440px the routes come before the plate. Primary
action: the two route links.

SIGNATURE: The place mark. The red number cell follows the reader through the
rail and the on-page contents. Plate 1 draws its trailing vortices aft once on
load. Motion grammar: state changes only, 150 to 200 ms ease-out; nothing else
moves, and reduced motion removes the plate draw.

FORM: standards manual (1970s federal graphics standards manuals), candidate 6
of 7 on the ordered grounded list; seed:c0393be6147afb34.

FINISH: unreviewed and undocumented is unfinished; this build ends with the finish review, the verdict, DESIGN.md, and every shipping raster carrying its provenance

## Open decisions

None blocking. Dark mode ships as a variant of the same sheet; light is the default.
