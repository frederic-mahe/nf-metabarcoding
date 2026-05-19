#!/usr/bin/env python3
"""Merge sub-/super-string OTUs onto their masters.

For each pair found by `search_for_terminal_gaps` (vsearch
``--cluster_smallmem --id 1.0``), the lower-abundance pupil OTU is
folded into the higher-abundance master:

* sample counts: ``pupil[s]`` added to ``master[s]``
* ``total``:    ``pupil_total`` added to ``master_total``
* ``cloud``:    ``master_cloud += pupil_cloud + 1`` (the legacy
                "+1 per merged pupil" quirk, locked in by [S39]'s
                golden-file test)
* ``spread``:   recomputed from the count of non-zero merged
                sample columns

OTUs that are neither master nor pupil pass through unchanged.
The output preserves the input column layout; pass-through rows
appear in their original order, masters appear last in master-row
insertion order. The bash wrapper sorts the result by OTU number.

If any OTU appears in both the master set and the pupil set the
script aborts non-zero with a WARNING on stderr — that overlap
would otherwise silently lose data.

Characterization tests at
``tests/python/test_merge_sub_superstring_otus.py`` pin byte-exact
output. Refactors must keep them green.
"""

from __future__ import annotations

__author__ = "Frédéric Mahé <frederic.mahe@cirad.fr>"
__date__ = "2019/09/20"
__version__ = "$Revision: 2.0"

import argparse
import csv
import sys
from typing import Optional


# The OTU table's per-cluster metadata columns. Anything else in the
# header is a sample column.
METADATA_COLUMNS: frozenset[str] = frozenset({
    "OTU", "total", "cloud", "amplicon", "length",
    "abundance", "chimera", "spread", "quality", "sequence",
    "identity", "taxonomy", "references",
})


def parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    """Parse arguments from the command line."""
    parser = argparse.ArgumentParser(
        description=(
            "Merge sub-/super-string OTUs onto their masters using a "
            "vsearch self-cluster .uc file's H lines."
        ),
    )
    parser.add_argument(
        "-t", "--table",
        required=True,
        help="input occurrence table (TSV with a metadata + sample header)",
    )
    parser.add_argument(
        "-m", "--matches",
        required=True,
        help=(
            "OTU connexions (vsearch .uc H lines, or any TSV ending in "
            "<pupil>\\t<master>)"
        ),
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="output table path (pass /dev/stdout to write to stdout)",
    )
    return parser.parse_args(argv)


def parse_connexions(path: str) -> dict[str, str]:
    """Return ``{pupil: master}`` from the match file.

    Each row's last two tab-separated fields are taken as
    ``(pupil, master)`` — matches vsearch's ``.uc`` H-line layout
    (last two cols = query, target).
    """
    connexions: dict[str, str] = {}
    with open(path) as f:
        for line in f:
            pupil, master = line.rstrip("\n").split("\t")[-2:]
            connexions[pupil] = master
    return connexions


class OverlapError(ValueError):
    """Raised when an OTU appears as both master and pupil."""


def check_no_overlap(connexions: dict[str, str]) -> None:
    """Abort if any OTU is both a master and a pupil.

    Mirrors the legacy behaviour: prints WARNING + the offending
    set to stderr, then ``sys.exit(-1)``.
    """
    pupils = set(connexions.keys())
    masters = set(connexions.values())
    overlap = masters & pupils
    if overlap:
        print(
            "WARNING: there are OTUs common to hit and query columns",
            file=sys.stderr,
        )
        print(overlap, file=sys.stderr)
        sys.exit(-1)


def merge_otus(
    table_path: str,
    output_path: str,
    connexions: dict[str, str],
) -> None:
    """Stream the OTU table, merge pupils onto masters, write output."""
    pupils = set(connexions.keys())
    masters = set(connexions.values())

    with open(table_path) as src, open(output_path, "w") as dst:
        reader = csv.DictReader(src, delimiter="\t")
        if reader.fieldnames is None:
            return
        sample_names = set(reader.fieldnames) - METADATA_COLUMNS
        master_rows: dict[str, dict[str, object]] = {}
        writer = csv.DictWriter(
            dst, fieldnames=reader.fieldnames, delimiter="\t",
        )
        writer.writeheader()
        for row in reader:
            otu = row["OTU"]
            if otu not in masters and otu not in pupils:
                writer.writerow(row)
                continue
            if otu in masters:
                master_rows[otu] = dict(row)
                for sample in sample_names:
                    master_rows[otu][sample] = int(row[sample])
                master_rows[otu]["spread"] = int(row["spread"])
                master_rows[otu]["total"] = int(row["total"])
                master_rows[otu]["cloud"] = int(row["cloud"])
            if otu in pupils:
                master = connexions[otu]
                merged = master_rows[master]
                for sample in sample_names:
                    merged[sample] = int(merged[sample]) + int(row[sample])
                merged["spread"] = sum(
                    1 for s in sample_names if int(merged[s]) > 0
                )
                merged["total"] = int(merged["total"]) + int(row["total"])
                merged["cloud"] = int(merged["cloud"]) + int(row["cloud"]) + 1
        for merged_row in master_rows.values():
            writer.writerow(merged_row)


def main(argv: Optional[list[str]] = None) -> int:
    args = parse_args(argv)
    connexions = parse_connexions(args.matches)
    check_no_overlap(connexions)
    merge_otus(args.table, args.output, connexions)
    return 0


if __name__ == "__main__":
    sys.exit(main())
