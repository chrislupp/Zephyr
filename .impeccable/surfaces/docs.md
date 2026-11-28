---
version: 1
slug: "docs"
primary_target: "docs"
related_targets: []
---

# Surface brief: documentation site (`docs/`)

Scope: the whole Docusaurus site. The site root is the start page; every doc
page, the API book, and the Copyright and License page share its look.
Visitor mode: Persuade on the start page (`docs/src/pages/index.tsx`); Read on
every doc page.

## Reader and job

Aerospace engineers and researchers. On the start page they decide in seconds
whether Zephyr belongs in their work and take Get started or the install line.
On a doc page they run a first analysis or wire Zephyr into an MDO loop.

The maintainer rejected the Approach Plate build: too much on the start page,
far too focused on the single optimization case, not modern. Wanted: the craft
of simulation vendor sites (Flow360, Simcenter) and open science code sites
(SU2, PyVista), a short start page, and doc pages with a new look but the
standard structure (sidebar, page, on-page contents). Ruled out: retro or
nostalgic, a stock docs template, a SaaS marketing tone. Chosen in a safer-
register round: Post-processor, dark by default.

Constraints: doc page text stays as written; released snapshots under
`docs/versioned_docs/` are untouched; Copyright and License stays reachable
from every page; no crate version in prose; no invented benchmarks or
adoption claims; every image is the library's own output.

## Direction contract

THESIS: The start page is a post-processor's result window. Zephyr's own
solution on the uCRM transport wing fills a dark viewport and the offer sits
over it. It refuses the long scroll and the feature-card grid.

OWN-WORLD: Graphite ground (#0a0e13 to #16202b) in a soft radial falloff,
near-white ink, and the viridis colour map as the only saturated colour: the
data's own colour bar, from which the one accent (a viridis teal-green) is
taken for links and the current place. Host Grotesk for text and headings, Red Hat
Mono for code and numeric readouts. Small radii (8 to 10px), 1px cool borders, no
glow. Light is the variant: the same structure on a cool white ground,
except the start page's result window and its figure tiles, which stay dark in
both themes, as a post-processor's viewport does (decided by the build: the
renders are drawn for a dark ground and the window is the thesis).

STORY: The visitor sees a transport wing solved by Zephyr with its loading,
reads one line saying what the library is, and takes Get started or copies
`cargo add zephyr-vlm`. Below, three capabilities (analysis, gradients for
MDO, validation), each shown by real output, and one short code sample.

FIRST VIEWPORT: Transparent navbar over the viewport (name with a viridis
mark, Docs, API, Examples, version, GitHub). Headline at about 64px over two
lines, one-sentence lede, Get started (filled, light) and the install command
(outlined, mono, copies on click). The uCRM wing render fills the right and
centre, its trailing legs fading aft. A caption top right names the case.
Along the bottom, readouts (C_L, C_Di, checked against OptVL 2.1.1) and the
colour bar for panel loading.

SIGNATURE: The result window. On load the render resolves from the panel mesh
into the coloured solution once (wireframe fades to loading, under a second);
reduced motion shows the solution at once. Hover and focus states are quick
(150ms) colour and border changes.

FORM: Post-processor (simulation vendor and post-processing render windows),
safer-register round 1, grounded candidate 1 of 2; seed:dfc5d3e5ebca392a.

FINISH: unreviewed and undocumented is unfinished; this build ends with the finish review, the verdict, DESIGN.md, and every shipping raster carrying its provenance

## Open decisions

None blocking.
