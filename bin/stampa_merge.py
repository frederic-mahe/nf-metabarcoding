#!/usr/bin/env python3
"""Parse vsearch ``--userout`` hits and emit a stampa-style results TSV.

This is the single-file successor to the legacy
``tmp/stampa/stampa_merge.py``: the legacy script consumed a directory
of ``hits.*`` shards (one per slurm array task) and produced one
``results.*`` per shard. Nextflow runs Part C as a single vsearch
invocation (no array split — `[S49]`), so this port collapses the
loop into one input → one output stream.

Input format (vsearch ``--userout`` with
``--userfields query+id<iddef>+target``)::

    amplicon_<abundance>\\tidentity\\thit

`hit` is either ``*`` (no hit) or ``accession<SPACE>taxonomy`` where
``taxonomy`` is a ``|``-separated lineage.

Output format (`[S33]` / `[S49]` shape)::

    amplicon\\tabundance\\tidentity\\ttaxonomy\\treferences

The last-common-ancestor is computed across the top hits for each
amplicon by intersecting taxonomy levels: levels that match across
every top hit are kept verbatim; levels that disagree are rewritten
to ``*``.
"""

from __future__ import annotations

__author__ = "Frédéric Mahé <frederic.mahe@cirad.fr>"
__date__ = "2026/05/19"
__version__ = "$Revision: 3.0"

import argparse
import sys
from typing import Iterator, Optional


def parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Compute the last-common-ancestor across vsearch "
            "--usearch_global top hits and emit a stampa-style "
            "results TSV."
        ),
    )
    parser.add_argument(
        "hits",
        help="vsearch --userout hits file (query+id+target columns)",
    )
    return parser.parse_args(argv)


def last_common_ancestor(taxonomies: list[list[str]]) -> list[str]:
    """Compute the last common ancestor across a list of taxonomies."""
    if not taxonomies:
        return []
    if len(taxonomies) == 1:
        return taxonomies[0]
    lca: list[str] = []
    for level in zip(*taxonomies):
        lca.append(level[0] if len(set(level)) == 1 else "*")
    return lca


def parse_hits(path: str) -> Iterator[tuple[str, str, str, list[str], str]]:
    """Yield ``(amplicon, abundance, identity, taxonomy, accession)`` rows.

    ``taxonomy`` is the ``|``-split lineage; ``accession`` is the
    reference accession (or ``"No_hit"``).
    """
    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            amplicon_field, identity, hit = line.split("\t")
            amplicon, abundance = amplicon_field.rsplit("_", 1)
            if hit == "*":
                accession = "No_hit"
                taxonomy = ["No_hit"]
            else:
                try:
                    accession, taxonomy_str = hit.split(" ", 1)
                except ValueError:
                    print(
                        f"stampa_merge: malformed hit "
                        f"{amplicon!r} {identity!r} {hit!r}",
                        file=sys.stderr,
                    )
                    raise
                taxonomy = taxonomy_str.split("|")
            yield amplicon, abundance, identity, taxonomy, accession


def emit_results(path: str, out) -> None:
    """Group top hits by amplicon and emit the LCA-merged TSV rows."""
    current_amplicon: Optional[str] = None
    current_abundance: str = ""
    current_identity: str = ""
    taxonomies: list[list[str]] = []
    accessions: list[str] = []

    def flush() -> None:
        if current_amplicon is None:
            return
        lca = last_common_ancestor(taxonomies)
        print(
            current_amplicon,
            current_abundance,
            current_identity,
            "|".join(lca),
            ",".join(accessions),
            sep="\t",
            file=out,
        )

    for amplicon, abundance, identity, taxonomy, accession in parse_hits(path):
        if current_amplicon != amplicon:
            flush()
            current_amplicon = amplicon
            current_abundance = abundance
            current_identity = identity
            taxonomies = [taxonomy]
            accessions = [accession]
        else:
            taxonomies.append(taxonomy)
            accessions.append(accession)
    flush()


def main(argv: Optional[list[str]] = None) -> int:
    args = parse_args(argv)
    emit_results(args.hits, sys.stdout)
    return 0


if __name__ == "__main__":
    sys.exit(main())
