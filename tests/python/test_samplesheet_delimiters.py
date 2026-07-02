"""Reject samplesheet cells carrying a TSV delimiter ([S95]).

``bin/parse_samplesheet.py`` emits its normalized rows as TSV (columns
joined by TAB, rows by newline). A cell that itself contains a TAB,
newline, or carriage return — reachable through CSV quoting, e.g. a
quoted ``"a<TAB>b"`` path cell — would shift columns or inject a phantom
row into that output. Such a cell must abort at startup, naming the row
and column. The ``sample`` column is already constrained to the stricter
``[S93]`` charset; this closes the same hole in the path / ``run``
columns.
"""

# COVERAGE: [S95]

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

from parse_samplesheet import SamplesheetError, read_samplesheet


def _write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8", newline="")
    return path


# ---------- read_samplesheet rejects embedded delimiters ------------------

def test_embedded_tab_in_path_cell_rejected(tmp_path: Path) -> None:
    # A quoted cell keeps its embedded TAB through the CSV reader.
    sheet = _write(
        tmp_path / "sheet.csv",
        'sample,fastq_1\nA,"reads\t1.fastq.gz"\n',
    )
    with pytest.raises(SamplesheetError) as excinfo:
        read_samplesheet(sheet)
    msg = str(excinfo.value)
    assert "fastq_1" in msg
    assert "2" in msg  # the offending row


def test_embedded_newline_in_path_cell_rejected(tmp_path: Path) -> None:
    sheet = _write(
        tmp_path / "sheet.csv",
        'sample,fastq_1\nA,"reads\n1.fastq.gz"\n',
    )
    with pytest.raises(SamplesheetError) as excinfo:
        read_samplesheet(sheet)
    assert "fastq_1" in str(excinfo.value)


def test_embedded_cr_in_run_cell_rejected(tmp_path: Path) -> None:
    sheet = _write(
        tmp_path / "sheet.csv",
        'sample,fastq_1,run\nA,reads_1.fastq.gz,"run\r17"\n',
    )
    with pytest.raises(SamplesheetError) as excinfo:
        read_samplesheet(sheet)
    assert "run" in str(excinfo.value)


def test_clean_cells_are_accepted(tmp_path: Path) -> None:
    sheet = _write(
        tmp_path / "sheet.csv",
        "sample,fastq_1,fastq_2\nA,a_1.fastq.gz,a_2.fastq.gz\n",
    )
    fieldnames, rows = read_samplesheet(sheet)
    assert fieldnames == ["sample", "fastq_1", "fastq_2"]
    assert rows[0][1]["fastq_1"] == "a_1.fastq.gz"


# ---------- CLI aborts on an embedded delimiter ---------------------------

def test_cli_rejects_embedded_tab(bin_dir: Path, tmp_path: Path) -> None:
    sheet = _write(
        tmp_path / "sheet.csv",
        'sample,fastq_1\nA,"reads\t1.fastq.gz"\n',
    )
    result = subprocess.run(
        [sys.executable, str(bin_dir / "parse_samplesheet.py"), str(sheet)],
        capture_output=True,
        text=True,
    )
    assert result.returncode != 0
    assert "fastq_1" in result.stderr
