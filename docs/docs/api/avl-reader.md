---
sidebar_position: 10
---

# AvlReader

```cpp
#include "vlm/io/avl_reader.h"
```

Parser for Athena Vortex Lattice (AVL) geometry files.

## Methods

```cpp
Aircraft ReadFile(const std::string& filename);
Aircraft ParseString(const std::string& content,
                     const std::string& base_path = "");
void SetVerbose(bool verbose);
```

### ReadFile

Reads an `.avl` file and returns a fully configured `Aircraft` object with all surfaces, sections, control surfaces, and reference quantities parsed from the file.

### ParseString

Parses AVL-format content from a string. The optional `base_path` is used to resolve relative paths for airfoil data files.

## Supported AVL Features

- Surface definitions with paneling parameters
- Section definitions with chord, twist, and position
- Control surface definitions (CONTROL keyword)
- Y-duplication (YDUPLICATE keyword)
- Angle offsets and translations
- Airfoil coordinate files (AFILE keyword)
- Body definitions
- Reference quantities (Sref, Cref, Bref)

## Example

```cpp
AvlReader reader;
reader.SetVerbose(true);  // print parsing details

Aircraft aircraft = reader.ReadFile("transport.avl");

std::cout << "Title: " << aircraft.GetTitle() << std::endl;
std::cout << "Surfaces: " << aircraft.GetNumSurfaces() << std::endl;
std::cout << "Total panels: " << aircraft.ComputeTotalVortexCount() << std::endl;
```

## Mass File Reader

The `MassReader` class parses AVL mass property files:

```cpp
#include "vlm/io/mass_reader.h"

MassReader mass_reader;
MassData mass_data = mass_reader.ReadFile("aircraft.mass");
```
