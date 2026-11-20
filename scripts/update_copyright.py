#!/usr/bin/env python3
# Project Zephyr
#
# Copyright 2024-2026, Christopher A. Lupp
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This work has been cleared for public release, distribution unlimited, case
# number: AFRL-2026-2198.
#
# The views expressed are those of the authors and do not reflect the
# official guidance or position of the United States Government, the
# Department of Defense or of the United States Air Force.
#
# Statement from DoD: The Appearance of external hyperlinks does not
# constitute endorsement by the United States Department of Defense (DoD) of
# the linked websites, of the information, products, or services contained
# therein. The DoD does not exercise any editorial, security, or other
# control over the information you may find at these locations.

"""Update copyright year ranges and the Cargo package version."""

import argparse
import re
import sys
from datetime import datetime
from pathlib import Path

COPYRIGHT_PATTERN = re.compile(
    r"(Copyright\s+)(\d{4})(-\d{4})?(,\s*Christopher\s+A\.\s+Lupp)"
)
VERSION_PATTERN = re.compile(r'(?m)^(version\s*=\s*")[^"]+(")$')
SOURCE_EXTENSIONS = {".rs", ".py", ".md", ".toml"}
SKIP_DIRS = {"build", "target", "docs", "node_modules"}


def update_file(path: Path, year: str, version: str | None, dry_run: bool) -> bool:
    """Update copyright years and version in a single file. Returns True if modified."""
    text = path.read_text()

    def replace_match(m):
        start_year = m.group(2)
        suffix = m.group(4)
        if start_year == year:
            return f"{m.group(1)}{start_year}{suffix}"
        return f"{m.group(1)}{start_year}-{year}{suffix}"

    new_text = COPYRIGHT_PATTERN.sub(replace_match, text)

    if version and path.name == "Cargo.toml":
        new_text = VERSION_PATTERN.sub(rf"\g<1>{version}\g<2>", new_text, count=1)

    if new_text != text:
        if not dry_run:
            path.write_text(new_text)
        return True
    return False


def main():
    parser = argparse.ArgumentParser(
        description="Update copyright year ranges and version numbers."
    )
    parser.add_argument(
        "--year",
        default=str(datetime.now().year),
        help="Target end year (default: current year)",
    )
    parser.add_argument(
        "--version",
        default=None,
        help="Target version string (e.g., 1.2.3). If omitted, versions are not updated.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be changed without modifying files",
    )
    parser.add_argument(
        "root",
        nargs="?",
        default=".",
        help="Root directory to scan (default: current directory)",
    )
    args = parser.parse_args()

    root = Path(args.root).resolve()
    updated = []
    scanned = 0

    for path in sorted(root.rglob("*")):
        if any(part in SKIP_DIRS for part in path.parts):
            continue
        if path.suffix not in SOURCE_EXTENSIONS:
            continue
        scanned += 1
        if update_file(path, args.year, args.version, args.dry_run):
            updated.append(path.relative_to(root))

    if updated:
        action = "Would update" if args.dry_run else "Updated"
        for f in updated:
            print(f"  {action}: {f}")
    else:
        print("All files up to date.")

    print(f"\nScanned {scanned} files, {len(updated)} {'would be ' if args.dry_run else ''}updated.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
