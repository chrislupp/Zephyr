---
sidebar_position: 10
---

# AVL and Mass Readers

`AvlReader::read_path`/`read_str` parse AVL geometry, modifiers, airfoils,
bodies, controls, and design variables. Relative airfoil files resolve from
the geometry path. `MassReader` parses AVL units, constants, scale/offset
directives, components, and products of inertia.

Fallible reads return `PegasusError` with path and line context.
