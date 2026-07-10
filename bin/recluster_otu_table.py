#!/usr/bin/env python3
"""Fold re-clustered member OTUs onto their centroids (post-mumu).

Optional, terminal coarse-clustering pass for divergent markers
([S103]/[S104], D20). ``recluster_search`` runs
``vsearch --cluster_size`` on the post-mumu FASTA and emits the ``^H``
lines of the ``.uc`` stream; each H line's last two tab-separated
fields are ``<member>\\t<centroid>`` labels of the shape
``<amplicon>;size=N;``. For every such pair the lower-abundance member
OTU is folded into its centroid:

* sample counts: ``member[s]`` added to ``centroid[s]``
* ``total``:      ``member_total`` added to ``centroid_total``
* ``spread``:     recomputed from the count of non-zero merged
                  sample columns
* ``cloud``:      left at the post-mumu ``"NA"`` — it is already a
                  string, so (unlike ``merge_substring_otus.py``'s
                  ``+1`` quirk) there is no arithmetic here
* ``chimera``:    left at ``"N"``
* metadata (``amplicon``/``length``/``abundance``/``quality``/
  ``sequence``/``identity``/``taxonomy``/``references``): taken from
  the centroid, i.e. the most abundant member under
  ``--cluster_size --sizein`` (D-c)

Surviving centroid rows (every OTU that is not a member) are emitted in
their input order and renumbered ``1..N`` (D-d), so no downstream sort
is needed. The total read count is asserted conserved.

This is a cousin of ``bin/merge_substring_otus.py`` kept separate so
[S39]'s golden output stays byte-exact. If any OTU appears in both the
centroid set and the member set the script aborts non-zero with a
WARNING on stderr — that overlap would otherwise silently lose data
(a ``--cluster_size`` centroid is never itself a member, so the guard
is trivially satisfied in practice).

Characterization tests at
``tests/python/test_recluster_otu_table.py`` pin byte-exact output.
Refactors must keep them green.
"""

from __future__ import annotations

__author__ = "Frédéric Mahé <frederic.mahe@cirad.fr>"
__date__ = "2026/07/10"
__version__ = "$Revision: 1.0"

import argparse
import csv
import sys
from typing import Optional


# The OTU table's per-cluster metadata columns. Anything else in the
# header is a sample column. Matches merge_substring_otus.py's set.
METADATA_COLUMNS: frozenset[str] = frozenset({
    "OTU", "total", "cloud", "amplicon", "length",
    "abundance", "chimera", "spread", "quality", "sequence",
    "identity", "taxonomy", "references",
})


def parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    """Parse arguments from the command line."""
    parser = argparse.ArgumentParser(
        description=(
            "Fold re-clustered member OTUs onto their centroids using a "
            "vsearch --cluster_size .uc file's H lines."
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
            "<member>\\t<centroid>); labels may carry a ;size=N; suffix"
        ),
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="output table path (pass /dev/stdout to write to stdout)",
    )
    return parser.parse_args(argv)


def strip_size(label: str) -> str:
    """Drop a trailing ``;size=N;`` annotation from a vsearch label.

    ``vsearch --sizeout`` writes ``<amplicon>;size=N;`` labels into the
    ``.uc`` stream; the occurrence table keys on the bare amplicon
    (column ``amplicon``), so the annotation is removed before matching.
    A label without the annotation is returned unchanged.
    """
    index = label.find(";size=")
    return label[:index] if index != -1 else label


def parse_connexions(path: str) -> dict[str, str]:
    """Return ``{member: centroid}`` from the match file.

    Each row's last two tab-separated fields are taken as
    ``(member, centroid)`` — matches vsearch's ``.uc`` H-line layout
    (last two cols = query, target) — with the ``;size=N;`` annotation
    stripped from both.
    """
    connexions: dict[str, str] = {}
    with open(path) as f:
        for line in f:
            member, centroid = line.rstrip("\n").split("\t")[-2:]
            connexions[strip_size(member)] = strip_size(centroid)
    return connexions


class OverlapError(ValueError):
    """Raised when an OTU appears as both centroid and member."""


def check_no_overlap(connexions: dict[str, str]) -> None:
    """Abort if any OTU is both a centroid and a member.

    Mirrors ``merge_substring_otus.py``: prints WARNING + the offending
    set to stderr, then ``sys.exit(-1)``.
    """
    members = set(connexions.keys())
    centroids = set(connexions.values())
    overlap = centroids & members
    if overlap:
        print(
            "WARNING: there are OTUs common to centroid and member columns",
            file=sys.stderr,
        )
        print(overlap, file=sys.stderr)
        sys.exit(-1)


def recluster_otus(
    table_path: str,
    output_path: str,
    connexions: dict[str, str],
) -> None:
    """Stream the OTU table, fold members onto centroids, write output."""
    members = set(connexions.keys())

    with open(table_path) as src:
        reader = csv.DictReader(src, delimiter="\t")
        fieldnames = reader.fieldnames
        if fieldnames is None:
            open(output_path, "w").close()
            return
        sample_names = set(fieldnames) - METADATA_COLUMNS
        rows = list(reader)

    by_amplicon: dict[str, dict[str, object]] = {
        row["amplicon"]: row for row in rows
    }

    # Surviving rows are every non-member, in input order (D-d).
    surviving = [row for row in rows if row["amplicon"] not in members]

    # Convert the accumulating columns to int on the surviving rows.
    for row in surviving:
        for sample in sample_names:
            row[sample] = int(row[sample])
        row["total"] = int(row["total"])

    # Fold each member into its centroid (order-independent — a member
    # may sit above its centroid in the post-mumu table).
    reads_before = sum(int(row["total"]) for row in rows)
    for member in members:
        member_row = by_amplicon.get(member)
        if member_row is None:
            continue  # a .uc label with no table row — skip defensively
        centroid_row = by_amplicon[connexions[member]]
        for sample in sample_names:
            centroid_row[sample] = (
                int(centroid_row[sample]) + int(member_row[sample])
            )
        centroid_row["total"] = (
            int(centroid_row["total"]) + int(member_row["total"])
        )

    with open(output_path, "w") as dst:
        writer = csv.DictWriter(
            dst, fieldnames=fieldnames, delimiter="\t",
        )
        writer.writeheader()
        reads_after = 0
        for new_otu, row in enumerate(surviving, start=1):
            row["spread"] = sum(
                1 for s in sample_names if int(row[s]) > 0
            )
            row["OTU"] = new_otu
            reads_after += int(row["total"])
            writer.writerow(row)

    if reads_after != reads_before:
        print(
            "ERROR: read count changed during re-clustering "
            f"({reads_before} -> {reads_after})",
            file=sys.stderr,
        )
        sys.exit(1)


def main(argv: Optional[list[str]] = None) -> int:
    args = parse_args(argv)
    connexions = parse_connexions(args.matches)
    check_no_overlap(connexions)
    recluster_otus(args.table, args.output, connexions)
    return 0


if __name__ == "__main__":
    sys.exit(main())
