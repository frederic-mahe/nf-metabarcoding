"""Characterization tests for ``bin/build_part_b_read_counts.py``.

Golden-file gate for the per-sample Part B read/cluster tracking
summary ([S107]) — the Part B counterpart of Part A's ``[S86]``
summary. Because Part B is pooled, the per-sample view is reconstructed
from the ``.distr`` ([S29]) and the per-sample columns of the
intermediate OTU tables, not from the (project-wide) step logs.

The fixture under ``tests/data/part_b_read_counts/`` exercises:

* ``S1`` — loses no reads (``reads_kept == reads_in``) but its cluster
  count drops 3 → 2 → 1 across substring-merge ([S39]) and mumu ([S44]);
* ``S2`` — loses reads at the occurrence-table filter ([S35], its
  ``amp_c`` cluster is dropped: ``reads_kept 3 < reads_in 11``) while
  staying in one cluster throughout;
* ``S3`` — an empty sample ([S09]) present only in ``--samples``: an
  all-zero row;
* the ``Total`` row (column sums);
* the optional ``clusters_recluster`` column ([S105]): present with
  ``--recluster``, absent without it.
"""

# COVERAGE: [S107]

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

SCRIPT = (
    Path(__file__).resolve().parents[2]
    / "bin"
    / "build_part_b_read_counts.py"
)
FIXTURE = (
    Path(__file__).resolve().parents[1]
    / "data"
    / "part_b_read_counts"
)

BASE_ARGS = [
    "--samples", "S1,S2,S3",
    "--distr", str(FIXTURE / "distribution.distr"),
    "--filtered", str(FIXTURE / "filtered.table"),
    "--merged", str(FIXTURE / "merged.table"),
    "--mumu", str(FIXTURE / "mumu.table"),
]

EXPECTED_NO_RECLUSTER = (
    "samples\treads_in\treads_kept\t"
    "clusters_kept\tclusters_merged\tclusters_mumu\n"
    "S1\t17\t17\t3\t2\t1\n"
    "S2\t11\t3\t1\t1\t1\n"
    "S3\t0\t0\t0\t0\t0\n"
    "Total\t28\t20\t4\t3\t2\n"
)

EXPECTED_WITH_RECLUSTER = (
    "samples\treads_in\treads_kept\t"
    "clusters_kept\tclusters_merged\tclusters_mumu\tclusters_recluster\n"
    "S1\t17\t17\t3\t2\t1\t1\n"
    "S2\t11\t3\t1\t1\t1\t1\n"
    "S3\t0\t0\t0\t0\t0\t0\n"
    "Total\t28\t20\t4\t3\t2\t2\n"
)


def _run(extra: list[str]) -> str:
    result = subprocess.run(
        [sys.executable, str(SCRIPT), *BASE_ARGS, *extra],
        check=True,
        capture_output=True,
        text=True,
    )
    return result.stdout


def test_summary_without_recluster() -> None:
    assert _run([]) == EXPECTED_NO_RECLUSTER


def test_summary_with_recluster() -> None:
    extra = ["--recluster", str(FIXTURE / "recluster.table")]
    assert _run(extra) == EXPECTED_WITH_RECLUSTER


def test_output_to_file(tmp_path: Path) -> None:
    out = tmp_path / "proj_2_samples_read_counts.tsv"
    subprocess.run(
        [sys.executable, str(SCRIPT), *BASE_ARGS, "-o", str(out)],
        check=True,
    )
    assert out.read_text() == EXPECTED_NO_RECLUSTER
