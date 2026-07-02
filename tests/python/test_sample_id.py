"""Unit tests for the shared sample-ID charset validator ([S93]).

A sample ID becomes a shell token and an output-file basename in the
Part A / Part B process scripts (``!{sampleId}``). The shared
``bin/sample_id.py`` module rejects any ID outside the safe set
``[A-Za-z0-9._-]`` (and any ID starting with ``-`` or ``.``) so a
malicious or mistyped ID cannot inject a command or escape the work
directory. The rule is enforced identically on the two entry paths:
the ``--input`` samplesheet (``parse_samplesheet.py``) and folder
discovery (``discover_fastq.py`` / ``discover_fasta.py``).
"""

# COVERAGE: [S93]

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

from sample_id import (
    InvalidSampleIdError,
    is_valid_sample_id,
    validate_sample_id,
)


# ---------- pure-function check -------------------------------------------

VALID_IDS = [
    "A",
    "B_run17",
    "sample.1",
    "s-1",
    "_underscore_lead",
    "123",
    "a.b-c_d",
    "paired_merge_ok",
]

INVALID_IDS = [
    "bad;id",        # command separator
    "a b",           # whitespace
    "a\tb",          # tab (would corrupt the normalized TSV)
    "a/b",           # path separator
    "$(whoami)",     # command substitution
    "a`id`b",        # backtick substitution
    "a|b",           # pipe
    "a&b",           # background / chain
    "a>b",           # redirection
    "a*b",           # glob
    "-leading",      # reads as an option to downstream tools
    ".hidden",       # hidden file
    "..",            # parent-dir traversal
    "",              # empty
    "café",          # non-ASCII
]


@pytest.mark.parametrize("sample_id", VALID_IDS)
def test_valid_ids_pass(sample_id: str) -> None:
    assert is_valid_sample_id(sample_id) is True
    # validate_sample_id returns the value unchanged when valid.
    assert validate_sample_id(sample_id) == sample_id


@pytest.mark.parametrize("sample_id", INVALID_IDS)
def test_invalid_ids_rejected(sample_id: str) -> None:
    assert is_valid_sample_id(sample_id) is False
    with pytest.raises(InvalidSampleIdError):
        validate_sample_id(sample_id)


def test_error_names_the_offending_id() -> None:
    with pytest.raises(InvalidSampleIdError, match="bad;id"):
        validate_sample_id("bad;id")


def test_context_is_included_in_message() -> None:
    with pytest.raises(InvalidSampleIdError, match="row 4"):
        validate_sample_id("bad;id", context="row 4")


# ---------- CLI integration: parse_samplesheet.py -------------------------

def _run(script: Path, *args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(script), *args],
        capture_output=True,
        text=True,
    )


def test_parse_samplesheet_rejects_bad_sample_id(
    bin_dir: Path, tmp_path: Path
) -> None:
    sheet = tmp_path / "sheet.csv"
    sheet.write_text("sample,fastq_1\nbad;id,reads_1.fastq.gz\n")
    result = _run(bin_dir / "parse_samplesheet.py", str(sheet))
    assert result.returncode != 0
    assert "bad;id" in result.stderr


def test_parse_samplesheet_accepts_clean_sample_id(
    bin_dir: Path, tmp_path: Path
) -> None:
    sheet = tmp_path / "sheet.csv"
    sheet.write_text("sample,fastq_1\ngood_1,reads_1.fastq.gz\n")
    result = _run(bin_dir / "parse_samplesheet.py", str(sheet))
    assert result.returncode == 0
    assert "good_1" in result.stdout


# ---------- CLI integration: discover_fastq.py ----------------------------

def test_discover_fastq_rejects_bad_sample_id(
    bin_dir: Path, tmp_path: Path
) -> None:
    # A single-end file whose stripped basename is an unsafe ID.
    (tmp_path / "bad;id.fastq").touch()
    result = _run(bin_dir / "discover_fastq.py", str(tmp_path))
    assert result.returncode != 0
    assert "bad;id" in result.stderr


# ---------- CLI integration: discover_fasta.py ----------------------------

def test_discover_fasta_rejects_bad_sample_id(
    bin_dir: Path, tmp_path: Path
) -> None:
    (tmp_path / "bad;id.fas").touch()
    result = _run(bin_dir / "discover_fasta.py", str(tmp_path))
    assert result.returncode != 0
    assert "bad;id" in result.stderr
