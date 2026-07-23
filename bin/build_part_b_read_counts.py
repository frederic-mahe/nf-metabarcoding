#!/usr/bin/env python3
"""Per-sample read/cluster tracking summary for Part B ([S107]).

The Part B counterpart of Part A's ``build_read_counts.sh`` ([S86]).
Part B is a *pooled* pipeline — every sample is merged at
``global_dereplication`` ([S31]) — so its step logs ([S45]) carry only
project-wide totals and cannot yield a per-sample view. The per-sample
dimension survives instead in the abundance data, so this summary is
reconstructed from:

* the sequence-to-sample distribution (``.distr``, [S29]): summing the
  ``size`` column per sample gives ``reads_in`` — the reads entering
  Part B;
* the per-sample columns (14…N of the occurrence-table schema) of the
  intermediate OTU tables at each curation stage: summing a sample's
  column gives its surviving read count, counting the rows where that
  column is ``> 0`` gives the number of clusters the sample appears in.

Reads are dropped at exactly one Part B step — the occurrence-table
filter ([S35]); every later curation step ([S39]/[S44]/[S104]) conserves
the read count and only reduces the number of clusters a sample belongs
to. Hence a single ``reads_kept`` column (from the post-filter table)
plus a ``clusters_*`` column per curation stage.

Output is a tab-separated table, one row per sample (sorted) then a
``Total`` row, with these columns::

    samples  reads_in  reads_kept
    clusters_kept  clusters_merged  clusters_mumu  [clusters_recluster]

Empty samples ([S09]) contribute no ``.distr`` rows and may carry no
table column; the authoritative ``--samples`` list is threaded in so
they still surface as an all-zero row.

Characterization tests at
``tests/python/test_build_part_b_read_counts.py`` pin the output.
"""

from __future__ import annotations

__author__ = "Frédéric Mahé <frederic.mahe@cirad.fr>"
__date__ = "2026/07/23"
__version__ = "$Revision: 1.0"

import argparse
import csv
import sys
from typing import Optional, TextIO

# The OTU table's per-cluster metadata columns ([S107] refers to the
# occurrence-table schema). Anything else in the header is a per-sample
# column. Matches merge_substring_otus.py / recluster_otu_table.py.
METADATA_COLUMNS: frozenset[str] = frozenset({
    "OTU", "total", "cloud", "amplicon", "length",
    "abundance", "chimera", "spread", "quality", "sequence",
    "identity", "taxonomy", "references",
})


def parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    """Parse arguments from the command line."""
    parser = argparse.ArgumentParser(
        description=(
            "Build the per-sample Part B read/cluster tracking summary "
            "([S107]) from the .distr and the intermediate OTU tables."
        ),
    )
    parser.add_argument(
        "--samples",
        required=True,
        help=(
            "comma-separated authoritative sample IDs; empty samples "
            "([S09]) still get a zero row"
        ),
    )
    parser.add_argument(
        "--distr",
        required=True,
        help="sequence-to-sample distribution file ([S29])",
    )
    parser.add_argument(
        "--filtered",
        required=True,
        help="post-filter occurrence table (build_occurrence_table, [S35])",
    )
    parser.add_argument(
        "--merged",
        required=True,
        help="post-substring-merge table (merge_substring_otus, [S39])",
    )
    parser.add_argument(
        "--mumu",
        required=True,
        help="post-mumu table (rebuild_post_mumu_table, [S44])",
    )
    parser.add_argument(
        "--recluster",
        default=None,
        help=(
            "reclustered table (recluster_merge, [S104]); when given, a "
            "clusters_recluster column is appended"
        ),
    )
    parser.add_argument(
        "-o", "--output",
        default="-",
        help="output path (default: stdout)",
    )
    return parser.parse_args(argv)


def parse_samples(raw: str) -> list[str]:
    """Return the sorted, de-duplicated, non-empty sample IDs."""
    return sorted({token.strip() for token in raw.split(",") if token.strip()})


def distr_reads(path: str) -> dict[str, int]:
    """Return ``{sample: reads_in}`` from a ``.distr`` file ([S29]).

    Each row is ``<sha1>\\t<sampleId>\\t<size>``; the per-sample reads
    entering Part B is the sum of the ``size`` column.
    """
    reads: dict[str, int] = {}
    with open(path) as handle:
        for line in handle:
            stripped = line.rstrip("\n")
            if not stripped:
                continue
            fields = stripped.split("\t")
            sample = fields[1]
            reads[sample] = reads.get(sample, 0) + int(fields[2])
    return reads


def table_reads_and_clusters(
    path: str,
) -> tuple[dict[str, int], dict[str, int]]:
    """Return ``(reads, clusters)`` per sample for one OTU table.

    ``reads[sample]`` is the sum of the sample's per-sample column;
    ``clusters[sample]`` is the number of rows where that column is
    ``> 0`` (the clusters the sample appears in). A header-only table
    yields zeros for every sample column it declares.
    """
    reads: dict[str, int] = {}
    clusters: dict[str, int] = {}
    with open(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = reader.fieldnames
        if not fieldnames:
            return reads, clusters
        samples = [col for col in fieldnames if col not in METADATA_COLUMNS]
        for sample in samples:
            reads[sample] = 0
            clusters[sample] = 0
        for row in reader:
            for sample in samples:
                value = int(row[sample])
                if value > 0:
                    reads[sample] += value
                    clusters[sample] += 1
    return reads, clusters


def build_rows(
    samples: list[str],
    reads_in: dict[str, int],
    filtered: tuple[dict[str, int], dict[str, int]],
    merged: tuple[dict[str, int], dict[str, int]],
    mumu: tuple[dict[str, int], dict[str, int]],
    recluster: Optional[tuple[dict[str, int], dict[str, int]]],
) -> tuple[list[str], list[list[str]]]:
    """Assemble the header and the per-sample + Total rows."""
    filtered_reads, filtered_clusters = filtered
    _, merged_clusters = merged
    _, mumu_clusters = mumu

    header = [
        "samples", "reads_in", "reads_kept",
        "clusters_kept", "clusters_merged", "clusters_mumu",
    ]
    if recluster is not None:
        header.append("clusters_recluster")

    _, recluster_clusters = recluster if recluster is not None else ({}, {})

    rows: list[list[str]] = []
    totals = [0] * (len(header) - 1)
    for sample in samples:
        values = [
            reads_in.get(sample, 0),
            filtered_reads.get(sample, 0),
            filtered_clusters.get(sample, 0),
            merged_clusters.get(sample, 0),
            mumu_clusters.get(sample, 0),
        ]
        if recluster is not None:
            values.append(recluster_clusters.get(sample, 0))
        rows.append([sample] + [str(value) for value in values])
        totals = [running + value for running, value in zip(totals, values)]

    rows.append(["Total"] + [str(value) for value in totals])
    return header, rows


def write_table(
    header: list[str],
    rows: list[list[str]],
    handle: TextIO,
) -> None:
    """Write the summary as a tab-separated table."""
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(header)
    writer.writerows(rows)


def main(argv: Optional[list[str]] = None) -> int:
    args = parse_args(argv)
    samples = parse_samples(args.samples)
    reads_in = distr_reads(args.distr)
    filtered = table_reads_and_clusters(args.filtered)
    merged = table_reads_and_clusters(args.merged)
    mumu = table_reads_and_clusters(args.mumu)
    recluster = (
        table_reads_and_clusters(args.recluster)
        if args.recluster is not None
        else None
    )

    header, rows = build_rows(
        samples, reads_in, filtered, merged, mumu, recluster,
    )

    if args.output == "-":
        write_table(header, rows, sys.stdout)
    else:
        with open(args.output, "w") as handle:
            write_table(header, rows, handle)
    return 0


if __name__ == "__main__":
    sys.exit(main())
