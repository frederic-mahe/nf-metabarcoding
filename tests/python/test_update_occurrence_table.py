"""Tests for bin/update_occurrence_table.py.

The join step must tolerate an assignments TSV that does or does not
carry a header row: the stampa path's published
``_taxonomy_stampa.tsv`` carries one ([S49]); the sintax
``_assignments_sintax.tsv`` intermediate and the empty-input fallback
do not. A leading row whose first field is the literal ``amplicon``
column name is a header sentinel (real amplicon IDs are SHA1 hashes,
[S46]) and is skipped rather than parsed as data.
"""

# COVERAGE: [S51]

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

from update_occurrence_table import main, parse_assignments

SCRIPT = (
    Path(__file__).resolve().parents[2] / "bin" / "update_occurrence_table.py"
)

HEADER = "amplicon\tabundance\tidentity\ttaxonomy\treferences\n"
ROW = "amp1\t10\t99.0\tBacteria|Firmicutes\trefA\n"


def test_parse_assignments_without_header(tmp_path: Path) -> None:
    path = tmp_path / "headerless.tsv"
    path.write_text(ROW)
    assert parse_assignments(str(path)) == {
        "amp1": ("99.0", "Bacteria|Firmicutes", "refA")
    }


def test_parse_assignments_skips_header(tmp_path: Path) -> None:
    path = tmp_path / "headered.tsv"
    path.write_text(HEADER + ROW)
    # The header row must not become an "amplicon"-keyed entry, and the
    # genuine amp1 row must survive.
    assert parse_assignments(str(path)) == {
        "amp1": ("99.0", "Bacteria|Firmicutes", "refA")
    }


def test_parse_assignments_empty_file(tmp_path: Path) -> None:
    path = tmp_path / "empty.tsv"
    path.write_text("")
    assert parse_assignments(str(path)) == {}


def _table() -> str:
    return (
        "OTU\ttotal\tcloud\tamplicon\tlength\tabundance\tchimera\tspread\t"
        "quality\tsequence\tidentity\ttaxonomy\treferences\ts1\n"
        "1\t10\tNA\tamp1\t10\t10\tN\t1\t0.001\tACGTACGTAC\t0.0\tNA\tNA\t10\n"
    )


def test_cli_splices_with_a_headered_assignments_file(
    tmp_path: Path, capsys
) -> None:
    table = tmp_path / "table.tsv"
    table.write_text(_table())
    assignments = tmp_path / "assignments.tsv"
    assignments.write_text(HEADER + ROW)

    rc = main(
        ["--occurrence_table", str(table), "--assignments", str(assignments)]
    )
    assert rc == 0

    out_lines = capsys.readouterr().out.splitlines()
    assert len(out_lines) == 2  # header + 1 data row
    fields = out_lines[1].split("\t")
    assert fields[3] == "amp1"
    assert (fields[10], fields[11], fields[12]) == (
        "99.0",
        "Bacteria|Firmicutes",
        "refA",
    )


def test_cli_runs_as_a_subprocess(tmp_path: Path) -> None:
    table = tmp_path / "table.tsv"
    table.write_text(_table())
    assignments = tmp_path / "assignments.tsv"
    assignments.write_text(HEADER + ROW)

    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--occurrence_table",
            str(table),
            "--assignments",
            str(assignments),
        ],
        capture_output=True,
        text=True,
        check=True,
    )
    rows = result.stdout.splitlines()
    assert len(rows) == 2
    assert rows[1].split("\t")[11] == "Bacteria|Firmicutes"
