"""Tests for bin/majority_assignment.py.

The pure helpers are unit-tested directly (importable module); the
CLI is pinned by a golden-file characterization test that reproduces
the byte-exact output of the legacy ``majority_assignment.py`` on a
fixture exercising every branch, across plain / ``.gz`` / ``.bz2``
reference inputs.
"""

# COVERAGE: [S66]

from __future__ import annotations

import bz2
import gzip
import subprocess
import sys
from pathlib import Path

from majority_assignment import (
    find_most_frequent_path,
    majority_for_references,
    open_reference,
    parse_fasta,
)


SCRIPT = Path(__file__).resolve().parents[2] / "bin" / "majority_assignment.py"


# ---------------------------------------------------------------------------
# Shared fixture (mirrors the legacy golden capture)
# ---------------------------------------------------------------------------

REFERENCE = (
    ">AB001 Bacteria|Firmicutes|Bacilli|Lactobacillales\n"
    "ACGTACGTACGTACGTACGT\n"
    ">AB002 Bacteria|Firmicutes|Bacilli|Bacillales\n"
    "ACGTACGTACGTACGTACGA\n"
    ">AB003 Bacteria|Firmicutes|Clostridia\n"
    "ACGTACGTACGTACGTACGC\n"
    ">AB004 Bacteria|Proteobacteria|Gammaproteobacteria|Enterobacterales\n"
    "ACGTACGTACGTACGTACGG\n"
)

INPUT_TABLE = (
    "OTU\ttotal\tcloud\tamplicon\tlength\tabundance\tchimera\tspread\t"
    "quality\tsequence\tidentity\ttaxonomy\treferences\ts1\n"
    "1\t10\tNA\tamp1\t20\t10\tN\t1\t0.001\tACGTACGTACGTACGTACGT\t99.0\t"
    "X\tAB001\t10\n"
    "2\t8\tNA\tamp2\t20\t8\tN\t1\t0.001\tACGTACGTACGTACGTACGA\t98.0\t"
    "X\tAB001,AB002\t8\n"
    "3\t6\tNA\tamp3\t20\t6\tN\t1\t0.001\tACGTACGTACGTACGTACGC\t97.0\t"
    "X\tAB001,AB002,AB003\t6\n"
    "4\t4\tNA\tamp4\t20\t4\tN\t1\t0.001\tACGTACGTACGTACGTACGG\t96.0\t"
    "X\tAB001,AB004\t4\n"
    "5\t2\tNA\tamp5\t20\t2\tN\t1\t0.001\tACGTACGTACGTACGTACGT\t0.0\t"
    "NA\tNo_hit\t2\n"
)

# Byte-exact output of the legacy script on the fixture above.
GOLDEN = (
    "OTU\tamplicon\ttaxonomy_majority\n"
    "1\tamp1\tBacteria (1/1)|Firmicutes (1/1)|Bacilli (1/1)|"
    "Lactobacillales (1/1)\n"
    "2\tamp2\tBacteria (2/2)|Firmicutes (2/2)|Bacilli (2/2)|"
    "Lactobacillales (1/2)\n"
    "3\tamp3\tBacteria (3/3)|Firmicutes (3/3)|Bacilli (2/3)|"
    "Lactobacillales (1/3)\n"
    "4\tamp4\tBacteria (2/2)|Firmicutes (1/2)|Bacilli (1/2)|"
    "Lactobacillales (1/2)\n"
    "5\tamp5\tNo_hit\n"
)


# ---------------------------------------------------------------------------
# find_most_frequent_path — pure helper
# ---------------------------------------------------------------------------

def test_single_path_each_rank_is_one_over_one() -> None:
    assert find_most_frequent_path([["A", "B", "C"]]) == \
        "A (1/1)|B (1/1)|C (1/1)"


def test_agreeing_paths_full_support() -> None:
    paths = [["A", "B", "C"], ["A", "B", "C"]]
    assert find_most_frequent_path(paths) == "A (2/2)|B (2/2)|C (2/2)"


def test_majority_wins_at_disagreeing_rank() -> None:
    # rank 2: B, B, X -> B wins 2/3.
    paths = [["A", "B"], ["A", "B"], ["A", "X"]]
    assert find_most_frequent_path(paths) == "A (3/3)|B (2/3)"


def test_tie_broken_by_first_seen_order() -> None:
    # rank 1: B vs X, tie -> first-seen (B) wins.
    paths = [["A", "B"], ["A", "X"]]
    assert find_most_frequent_path(paths) == "A (2/2)|B (1/2)"


def test_unequal_depth_padded_with_na() -> None:
    # The shallower lineage contributes "NA" at the deepest rank;
    # at rank 3 the count is {NA: 1, D: 1}, tie -> first-seen "NA".
    paths = [["A", "B", "C"], ["A", "B", "C", "D"]]
    assert find_most_frequent_path(paths) == \
        "A (2/2)|B (2/2)|C (2/2)|NA (1/2)"


# ---------------------------------------------------------------------------
# majority_for_references
# ---------------------------------------------------------------------------

ACCESSIONS = {
    "AB001": "Bacteria|Firmicutes|Bacilli",
    "AB002": "Bacteria|Firmicutes|Clostridia",
}


def test_no_hit_passthrough() -> None:
    assert majority_for_references("No_hit", ACCESSIONS) == "No_hit"


def test_single_reference() -> None:
    assert majority_for_references("AB001", ACCESSIONS) == \
        "Bacteria (1/1)|Firmicutes (1/1)|Bacilli (1/1)"


def test_multiple_references_disagree() -> None:
    assert majority_for_references("AB001,AB002", ACCESSIONS) == \
        "Bacteria (2/2)|Firmicutes (2/2)|Bacilli (1/2)"


# ---------------------------------------------------------------------------
# open_reference / parse_fasta — compression handling (modernization)
# ---------------------------------------------------------------------------

def _write_plain(tmp_path: Path) -> Path:
    path = tmp_path / "ref.fas"
    path.write_text(REFERENCE)
    return path


def _write_gz(tmp_path: Path) -> Path:
    path = tmp_path / "ref.fas.gz"
    with gzip.open(path, "wt") as handle:
        handle.write(REFERENCE)
    return path


def _write_bz2(tmp_path: Path) -> Path:
    path = tmp_path / "ref.fas.bz2"
    with bz2.open(path, "wt") as handle:
        handle.write(REFERENCE)
    return path


def test_open_reference_reads_each_codec(tmp_path: Path) -> None:
    for builder in (_write_plain, _write_gz, _write_bz2):
        path = builder(tmp_path)
        with open_reference(str(path)) as handle:
            assert handle.read() == REFERENCE


def test_parse_fasta_indexes_accession_to_lineage(tmp_path: Path) -> None:
    accessions = parse_fasta(str(_write_gz(tmp_path)))
    assert accessions["AB001"] == "Bacteria|Firmicutes|Bacilli|Lactobacillales"
    assert accessions["AB003"] == "Bacteria|Firmicutes|Clostridia"
    # only header lines are indexed.
    assert len(accessions) == 4


def test_parse_fasta_codec_agnostic(tmp_path: Path) -> None:
    plain = parse_fasta(str(_write_plain(tmp_path)))
    gz = parse_fasta(str(_write_gz(tmp_path)))
    bz = parse_fasta(str(_write_bz2(tmp_path)))
    assert plain == gz == bz


# ---------------------------------------------------------------------------
# CLI golden-file characterization across the three reference codecs
# ---------------------------------------------------------------------------

def _run(input_table: Path, reference: Path) -> str:
    completed = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--input_table",
            str(input_table),
            "--reference_db",
            str(reference),
        ],
        capture_output=True,
        text=True,
        check=True,
    )
    return completed.stdout


def test_cli_matches_golden_plain(tmp_path: Path) -> None:
    table = tmp_path / "table.tsv"
    table.write_text(INPUT_TABLE)
    assert _run(table, _write_plain(tmp_path)) == GOLDEN


def test_cli_matches_golden_gz(tmp_path: Path) -> None:
    table = tmp_path / "table.tsv"
    table.write_text(INPUT_TABLE)
    assert _run(table, _write_gz(tmp_path)) == GOLDEN


def test_cli_matches_golden_bz2(tmp_path: Path) -> None:
    table = tmp_path / "table.tsv"
    table.write_text(INPUT_TABLE)
    assert _run(table, _write_bz2(tmp_path)) == GOLDEN
