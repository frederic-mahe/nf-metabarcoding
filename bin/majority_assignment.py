#!/usr/bin/env python3
"""Compute a majority-rule taxonomic assignment per OTU ([S66]).

This is the modernized port of the legacy ``majority_assignment.py``
(F. Mahé, 2024). It runs as the opt-in final step of the regular
Part C path: after ``update_occurrence_table`` ([S51]) has spliced the
stampa assignments onto the occurrence table, this script revisits the
per-OTU list of reference accessions (the ``references`` column) and,
for each taxonomic rank, reports the **most frequent** name across the
lineages of those accessions, annotated with its support.

Inputs:
    * ``--input_table`` — the regular assigned occurrence table
      (``<basename>_table_assigned.tsv``). Only the ``OTU``,
      ``amplicon``, and ``references`` columns are consumed; the table
      is read with :class:`csv.DictReader`, so column order is
      irrelevant as long as those three headers are present.
    * ``--reference_db`` — the stampa-formatted reference fasta used
      for the assignment ([S47]). Each header is
      ``>accession<SPACE>lineage`` where ``lineage`` is a
      ``|``-separated rank path. Plain, gzip- (``.gz``) and bzip2-
      (``.bz2``) compressed fasta are all accepted (detected by magic
      bytes, not by extension).

Output (written to stdout, three tab-separated columns)::

    OTU\tamplicon\ttaxonomy_majority

For each rank the majority name is formatted ``name (hits/total)``
where ``total`` is the number of reference accessions for the OTU and
``hits`` is how many of them carry ``name`` at that rank; ranks are
``|``-joined. Lineages of unequal depth are padded with ``NA``
(:func:`itertools.zip_longest`), so a shallower lineage contributes
``NA`` at the deeper ranks. Ties are broken by first-seen order
(:meth:`collections.Counter.most_common`). A row whose ``references``
is ``No_hit`` is passed through verbatim as ``No_hit``.

Majority assignment only makes sense for the stampa method: sintax
([S50]) leaves ``references`` at the ``NA`` placeholder, so the
workflow forbids the flag combination at startup ([S66]).
"""

from __future__ import annotations

__authors__ = "Frédéric Mahé"
__license__ = "GPL3"
__email__ = "frederic.mahe@cirad.fr"
__version__ = "0.3"

import bz2
import csv
import gzip
import itertools
import sys
from argparse import ArgumentParser, Namespace
from collections import Counter
from typing import IO, Iterator, Optional, TextIO


# gzip and bzip2 file signatures (first two bytes).
GZIP_MAGIC = b"\x1f\x8b"
BZIP2_MAGIC = b"BZ"


def parse_args(argv: Optional[list[str]] = None) -> Namespace:
    """Parse arguments from the command line."""
    parser = ArgumentParser(
        description="parse occurrence table and compute majority assignments."
    )
    parser.add_argument(
        "-i",
        "--input_table",
        action="store",
        nargs="?",
        dest="input_table",
        required=True,
        help="assigned occurrence table (<basename>_table_assigned.tsv)",
    )
    parser.add_argument(
        "-r",
        "--reference_db",
        action="store",
        nargs="?",
        dest="reference_db",
        required=True,
        help="stampa-formatted reference fasta (plain, .gz or .bz2)",
    )
    return parser.parse_args(argv)


def open_reference(filename: str) -> TextIO:
    """Open a fasta reference as text, transparently handling compression.

    The codec is chosen from the file's magic bytes rather than its
    extension so a mislabelled or extension-less file still opens
    correctly.
    """
    with open(filename, "rb") as probe:
        signature = probe.read(2)
    if signature == GZIP_MAGIC:
        return gzip.open(filename, "rt")
    if signature == BZIP2_MAGIC:
        return bz2.open(filename, "rt")
    return open(filename, "rt")


def parse_fasta(filename: str) -> dict[str, str]:
    """Read a reference fasta and index accession -> lineage string."""
    accessions: dict[str, str] = {}
    with open_reference(filename) as input_file:
        for line in input_file:
            if not line.startswith(">"):
                continue
            accession, taxonomy = line.strip().lstrip(">").split(" ", 1)
            accessions[accession] = taxonomy
    return accessions


def find_most_frequent_path(all_paths: list[list[str]]) -> str:
    """Report the most frequent name at each rank, with its support.

    ``all_paths`` is one ``|``-split lineage per reference accession.
    Ranks are walked with :func:`itertools.zip_longest` (fill value
    ``"NA"``) so unequal-depth lineages still align; the per-rank
    winner is ``Counter.most_common(1)`` (ties broken by first-seen
    order). The result is ``name (hits/total)`` per rank, ``|``-joined.
    """
    new_path: list[str] = []
    n_hits = len(all_paths)
    for taxo_level in itertools.zip_longest(*all_paths, fillvalue="NA"):
        count = Counter(taxo_level)
        most_common = count.most_common(1)[0]
        taxo_name, occurrences = most_common
        new_path.append(f"{taxo_name} ({occurrences}/{n_hits})")
    return "|".join(new_path)


def majority_for_references(
    references: str, accessions: dict[str, str]
) -> str:
    """Compute the majority lineage for one OTU's ``references`` field.

    Returns ``"No_hit"`` verbatim for an unassigned OTU; otherwise the
    ``|``-joined ``name (hits/total)`` path from
    :func:`find_most_frequent_path`.
    """
    if references == "No_hit":
        return "No_hit"
    assignments = [
        accessions[accession].split("|")
        for accession in references.split(",")
    ]
    return find_most_frequent_path(assignments)


def majority_rows(
    filename: str, accessions: dict[str, str]
) -> Iterator[tuple[str, str, str]]:
    """Yield ``(OTU, amplicon, taxonomy_majority)`` for each table row."""
    csv.field_size_limit(sys.maxsize)  # some reference lists are very long
    with open(filename) as csvfile:
        tablereader = csv.DictReader(csvfile, delimiter="\t")
        for row in tablereader:
            yield (
                row["OTU"],
                row["amplicon"],
                majority_for_references(row["references"], accessions),
            )


def emit_table(
    rows: Iterator[tuple[str, str, str]], out: IO[str] = sys.stdout
) -> None:
    """Write the majority-assignment table (header + rows) to ``out``."""
    print("OTU", "amplicon", "taxonomy_majority", sep="\t", file=out)
    for otu, amplicon, taxonomy_majority in rows:
        print(otu, amplicon, taxonomy_majority, sep="\t", file=out)


def main(argv: Optional[list[str]] = None) -> int:
    """Merge multiple top-hit assignments into a majority assignment."""
    args = parse_args(argv)
    accessions = parse_fasta(args.reference_db)
    emit_table(majority_rows(args.input_table, accessions))
    return 0


if __name__ == "__main__":
    sys.exit(main())
