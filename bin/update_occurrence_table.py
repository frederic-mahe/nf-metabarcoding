#!/usr/bin/env python3
"""Splice Part C's taxonomy assignments onto the occurrence table.

Skeleton implementation for `[S51]`: this is the join step that
overwrites the ``identity``, ``taxonomy``, and ``references`` columns
of Part B's `[S46]` occurrence table with the values produced by
``stampa_merge.py`` (or sintax, see `[S50]`).

The script is intentionally minimal — D04 still needs to resolve how
Part C publishes its output and which input modes are supported.
Until then, the script:

    * reads the occurrence table from ``--occurrence_table``;
    * reads the assignments TSV (one row per amplicon, columns
      ``amplicon\\tabundance\\tidentity\\ttaxonomy\\treferences``)
      from ``--assignments``;
    * writes the updated table to stdout, replacing columns 11
      (``identity``), 12 (``taxonomy``), and 13 (``references``).

Rows whose amplicon is missing from the assignments file are passed
through unchanged so an empty / partial assignments file does not
silently drop data.
"""

from __future__ import annotations

__author__ = "Frédéric Mahé <frederic.mahe@cirad.fr>"
__date__ = "2026/05/19"
__version__ = "$Revision: 0.1"

import argparse
import sys
from typing import Optional


# Column indices in the occurrence table ([S46]).
AMPLICON_COL = 3       # 0-based; column 4 in human terms
IDENTITY_COL = 10
TAXONOMY_COL = 11
REFERENCES_COL = 12


def parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Splice taxonomic assignments onto the Part B occurrence "
            "table by amplicon ID. Skeleton implementation — see "
            "DECISIONS.md D04."
        ),
    )
    parser.add_argument(
        "--occurrence_table",
        required=True,
        help="Part B occurrence table (<basename>_table.tsv)",
    )
    parser.add_argument(
        "--assignments",
        required=True,
        help="taxonomic assignments TSV "
             "(amplicon\\tabundance\\tidentity\\ttaxonomy\\treferences)",
    )
    return parser.parse_args(argv)


def parse_assignments(path: str) -> dict[str, tuple[str, str, str]]:
    """Map ``amplicon -> (identity, taxonomy, references)``."""
    assignments: dict[str, tuple[str, str, str]] = {}
    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            # amplicon\tabundance\tidentity\ttaxonomy\treferences
            amplicon, _abundance, identity, taxonomy, references = fields
            assignments[amplicon] = (identity, taxonomy, references)
    return assignments


def main(argv: Optional[list[str]] = None) -> int:
    args = parse_args(argv)
    assignments = parse_assignments(args.assignments)

    with open(args.occurrence_table) as handle:
        header = handle.readline().rstrip("\n")
        print(header)
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            amplicon = fields[AMPLICON_COL]
            update = assignments.get(amplicon)
            if update is not None:
                identity, taxonomy, references = update
                fields[IDENTITY_COL] = identity
                fields[TAXONOMY_COL] = taxonomy
                fields[REFERENCES_COL] = references
            print("\t".join(fields))
    return 0


if __name__ == "__main__":
    sys.exit(main())
