"""Unit tests for ``bin/parse_samplesheet.py`` ([S70]).

The helper does structural validation of the ``--input`` samplesheet in
two profiles (fastq / fasta), inferred from the columns present. File
existence is intentionally *not* checked here — Nextflow enforces it
downstream via ``file(checkIfExists: true)`` — so every case below runs
from strings, no filesystem needed.
"""

# COVERAGE: [S70]

from __future__ import annotations

import os
import subprocess
import sys

import pytest

from parse_samplesheet import (
    FastaSample,
    FastqSample,
    SamplesheetError,
    infer_profile,
    parse_rows,
    to_tsv,
)


def _rows(*dicts: dict) -> list[tuple[int, dict]]:
    """Pair row dicts with synthetic 1-based-after-header line numbers."""
    return [(i, d) for i, d in enumerate(dicts, start=2)]


# --- profile inference -------------------------------------------------

def test_infer_profile_fastq() -> None:
    assert infer_profile(["sample", "fastq_1", "fastq_2"]) == "fastq"


def test_infer_profile_fasta() -> None:
    assert infer_profile(["sample", "fasta", "qual"]) == "fasta"


def test_infer_profile_unresolvable() -> None:
    with pytest.raises(SamplesheetError, match="cannot infer"):
        infer_profile(["sample", "reads"])


# --- valid sheets ------------------------------------------------------

def test_valid_fastq_paired_and_single_end() -> None:
    fields = ["sample", "fastq_1", "fastq_2", "run"]
    profile, records = parse_rows(
        fields,
        _rows(
            {"sample": "A", "fastq_1": "A_1.fq", "fastq_2": "A_2.fq",
             "run": "r17"},
            {"sample": "B", "fastq_1": "B.fq", "fastq_2": "", "run": ""},
        ),
    )
    assert profile == "fastq"
    assert records[0] == FastqSample("A", "A_1.fq", "A_2.fq", "r17")
    assert records[0].single_end is False
    assert records[1].single_end is True


def test_valid_fasta_with_sibling_defaulting() -> None:
    fields = ["sample", "fasta"]
    profile, records = parse_rows(
        fields, _rows({"sample": "A", "fasta": "/d/A.fas"})
    )
    assert profile == "fasta"
    assert records[0] == FastaSample(
        "A", "/d/A.fas", "/d/A.qual", "/d/A.stats"
    )


def test_fasta_explicit_qual_stats_override_siblings() -> None:
    fields = ["sample", "fasta", "qual", "stats"]
    _, records = parse_rows(
        fields,
        _rows({"sample": "A", "fasta": "/d/A.fas",
               "qual": "/q/x.qual", "stats": "/s/x.stats"}),
    )
    assert records[0].qual == "/q/x.qual"
    assert records[0].stats == "/s/x.stats"


def test_fasta_notmerged_sample_is_shadow() -> None:
    fields = ["sample", "fasta"]
    _, records = parse_rows(
        fields, _rows({"sample": "X_notmerged", "fasta": "/d/X_notmerged.fas"})
    )
    assert records[0].shadow is True


# --- structural errors -------------------------------------------------

def test_missing_required_column() -> None:
    # fastq profile (has fastq_1) but no `sample` column.
    with pytest.raises(SamplesheetError, match="missing required column"):
        parse_rows(["fastq_1"], _rows({"fastq_1": "A.fq"}))


def test_unknown_column_rejected() -> None:
    fields = ["sample", "fastq_1", "bogus"]
    with pytest.raises(SamplesheetError, match="unknown column"):
        parse_rows(fields, _rows({"sample": "A", "fastq_1": "A.fq",
                                  "bogus": "1"}))


def test_empty_required_cell() -> None:
    with pytest.raises(SamplesheetError, match="required column 'fastq_1'"):
        parse_rows(["sample", "fastq_1"],
                   _rows({"sample": "A", "fastq_1": ""}))


def test_duplicate_sample_lists_both_rows() -> None:
    fields = ["sample", "fastq_1"]
    rows = _rows(
        {"sample": "A", "fastq_1": "A_1.fq"},
        {"sample": "A", "fastq_1": "A_2.fq"},
    )
    with pytest.raises(SamplesheetError) as excinfo:
        parse_rows(fields, rows)
    message = str(excinfo.value)
    assert "duplicate sample" in message
    assert "2" in message and "3" in message  # both row numbers named


def test_reserved_notmerged_rejected_in_fastq_profile() -> None:
    fields = ["sample", "fastq_1"]
    with pytest.raises(SamplesheetError, match="reserved suffix"):
        parse_rows(fields, _rows({"sample": "X_notmerged",
                                  "fastq_1": "x.fq"}))


# --- tilde expansion ([S60]) ------------------------------------------

def test_tilde_expansion_in_paths() -> None:
    fields = ["sample", "fastq_1", "fastq_2"]
    _, records = parse_rows(
        fields,
        _rows({"sample": "A", "fastq_1": "~/A_1.fq", "fastq_2": "~/A_2.fq"}),
    )
    home = os.path.expanduser("~")
    assert records[0].fastq_1 == f"{home}/A_1.fq"
    assert records[0].fastq_2 == f"{home}/A_2.fq"


# --- normalized TSV + CLI ---------------------------------------------

def test_to_tsv_round_trips_header_and_rows() -> None:
    fields = ["sample", "fastq_1", "fastq_2", "run"]
    profile, records = parse_rows(
        fields,
        _rows({"sample": "A", "fastq_1": "A_1.fq", "fastq_2": "", "run": ""}),
    )
    tsv = to_tsv(profile, records)
    lines = tsv.splitlines()
    assert lines[0] == "sample\tfastq_1\tfastq_2\trun"
    assert lines[1] == "A\tA_1.fq\t\t"


def test_cli_prints_tsv(tmp_path, bin_dir) -> None:
    sheet = tmp_path / "sheet.csv"
    sheet.write_text(
        "sample,fastq_1,fastq_2\n"
        "A,/d/A_1.fq,/d/A_2.fq\n"
        "B,/d/B.fq,\n"
    )
    proc = subprocess.run(
        [sys.executable, str(bin_dir / "parse_samplesheet.py"), str(sheet)],
        capture_output=True, text=True, check=True,
    )
    lines = proc.stdout.splitlines()
    assert lines[0] == "sample\tfastq_1\tfastq_2\trun"
    assert lines[1] == "A\t/d/A_1.fq\t/d/A_2.fq\t"
    assert lines[2] == "B\t/d/B.fq\t\t"


def test_cli_errors_nonzero_on_duplicate(tmp_path, bin_dir) -> None:
    sheet = tmp_path / "dupe.csv"
    sheet.write_text("sample,fasta\nA,/d/A.fas\nA,/d/A2.fas\n")
    proc = subprocess.run(
        [sys.executable, str(bin_dir / "parse_samplesheet.py"), str(sheet)],
        capture_output=True, text=True,
    )
    assert proc.returncode == 1
    assert "duplicate sample" in proc.stderr
