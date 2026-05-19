#!/usr/bin/env python3
"""Rebuild an occurrence table with 13 fields of metadata after mumu.

mumu (``--new_otu_table``) emits a stripped table — amplicon ID
plus one column per sample. This script re-attaches the per-amplicon
metadata (``length``, ``abundance``, ``quality``, ``sequence``,
``identity``, ``taxonomy``, ``references``) from the pre-mumu OTU
table, recomputes ``total`` (sum of mumu sample columns) and
``spread`` (count of non-zero mumu sample columns), replaces
``cloud`` with ``"NA"`` (mumu drops the cloud info), forces
``chimera`` to ``"N"``, and renumbers OTUs starting at 1.

The downstream ``size=0 → 1`` awk hotfix is applied by the
nextflow `rebuild_post_mumu_table` wrapper, not by this script.

Characterization tests at
``tests/python/test_rebuild_table_after_mumu.py`` pin byte-exact
output. Refactors must keep them green.
"""

from __future__ import annotations

__author__ = "Frédéric Mahé <frederic.mahe@cirad.fr>"
__date__ = "2020/12/22"
__version__ = "$Revision: 2.0"

import argparse
import sys
from dataclasses import dataclass
from typing import Optional


# Constants the legacy bash + python pair always uses.
CLOUD: str = "NA"        # cloud info is lost during the mumu pass
CHIMERA: str = "N"       # mumu only sees non-chimeric OTUs already


@dataclass(frozen=True)
class AmpliconMetadata:
    """Per-amplicon metadata sliced from the pre-mumu OTU table."""

    length: str
    abundance: str
    quality: str
    sequence: str
    identity: str
    taxonomy: str
    references: str


def parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    """Parse arguments from the command line."""
    parser = argparse.ArgumentParser(
        description=(
            "Rebuild an occurrence table with 13 fields of metadata "
            "from a mumu --new_otu_table output."
        ),
    )
    parser.add_argument(
        "--mumu_table",
        required=True,
        help="OTU table emitted by mumu (--new_otu_table)",
    )
    parser.add_argument(
        "--old_table",
        required=True,
        help="pre-mumu OTU table (carries the per-amplicon metadata)",
    )
    return parser.parse_args(argv)


def parse_old_table(
    path: str,
) -> tuple[str, dict[str, AmpliconMetadata]]:
    """Return ``(header_line, {amplicon: metadata})`` from the old table."""
    print("PROGRESS: parsing old table", file=sys.stderr)
    amplicons: dict[str, AmpliconMetadata] = {}
    with open(path) as f:
        header = next(f).rstrip("\n")
        for line in f:
            cols = line.rstrip("\n").split("\t")
            amplicons[cols[3]] = AmpliconMetadata(
                length=cols[4],
                abundance=cols[5],
                quality=cols[8],
                sequence=cols[9],
                identity=cols[10],
                taxonomy=cols[11],
                references=cols[12],
            )
    return header, amplicons


def write_rebuilt_rows(
    mumu_table_path: str,
    amplicons: dict[str, AmpliconMetadata],
) -> None:
    """Stream the mumu table and emit rebuilt rows on stdout."""
    print("PROGRESS: parsing mumu table", file=sys.stderr)
    with open(mumu_table_path) as f:
        next(f)  # discard mumu's header
        for line_counter, line in enumerate(f, start=1):
            cols = line.rstrip("\n").split("\t")
            amplicon = cols[0]
            samples = cols[1:]
            meta = amplicons[amplicon]
            print(
                line_counter,
                sum(int(c) for c in samples),
                CLOUD,
                amplicon,
                meta.length,
                meta.abundance,
                CHIMERA,
                sum(1 for c in samples if c != "0"),
                meta.quality,
                meta.sequence,
                meta.identity,
                meta.taxonomy,
                meta.references,
                "\t".join(samples),
                sep="\t",
            )


def main(argv: Optional[list[str]] = None) -> int:
    args = parse_args(argv)
    header, amplicons = parse_old_table(args.old_table)
    print(header)
    write_rebuilt_rows(args.mumu_table, amplicons)
    return 0


if __name__ == "__main__":
    sys.exit(main())
