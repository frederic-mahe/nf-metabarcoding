"""Tests for bin/check_reference.py, the full-reference check (D22)."""

# COVERAGE: [S115], [S116], [S117]

import bz2
import gzip
import subprocess
import sys
from pathlib import Path

import pytest

from check_reference import (
    ERRORS,
    WARNINGS,
    Finding,
    check_reference,
    status,
    summarise,
)

BIN = Path(__file__).resolve().parents[2] / "bin" / "check_reference.py"
SEQ = b"ACGTACGTACGTACGTACGTACGTACGTACGT"
STAMPA_OK = b">refA Bacteria|Proteobacteria|Gammaproteobacteria"
SINTAX_OK = b">refA;tax=d:Bacteria,p:Proteobacteria;"


def _reference(tmp_path: Path, headers: list[bytes],
               name: str = "reference.fas") -> Path:
    """Write a fasta whose records carry `headers`, one sequence each."""
    path = tmp_path / name
    path.write_bytes(b"".join(h + b"\n" + SEQ + b"\n" for h in headers))
    return path


def _issues(findings: list[Finding]) -> set[str]:
    return {f.issue for f in findings}


# ---------------------------------------------------------------------------
# issue codes
# ---------------------------------------------------------------------------

def test_issue_codes_split_into_errors_and_warnings() -> None:
    # [S115] / [S116]: each code is always an error or always a warning.
    assert ERRORS == {
        "format", "invalid_utf8", "control_char",
        "duplicate_accession", "comma_in_accession", "whitespace",
    }
    assert WARNINGS == {"empty_rank", "star_in_lineage"}


# ---------------------------------------------------------------------------
# clean references
# ---------------------------------------------------------------------------

@pytest.mark.parametrize(
    ("grammar", "headers"),
    [
        ("stampa", [STAMPA_OK, b">refB Bacteria|Firmicutes"]),
        ("sintax", [SINTAX_OK, b">refB;tax=d:Bacteria,p:Firmicutes"]),
        # hybrid: both annotations, read by both grammars ([S64])
        ("stampa", [b">refA;tax=d:Bacteria,p:Proteobacteria; "
                    b"Bacteria|Proteobacteria"]),
        ("sintax", [b">refA;tax=d:Bacteria,p:Proteobacteria; "
                    b"Bacteria|Proteobacteria"]),
        # non-ASCII valid UTF-8 and spaces inside a rank name are fine
        ("stampa", [STAMPA_OK,
                    ">refB Eukaryota|Homo sapiens|Genus_é".encode()]),
        # uneven depth across the file is not reported ([S116])
        ("stampa", [STAMPA_OK, b">refB Bacteria"]),
    ],
)
def test_clean_reference_reports_nothing(
    tmp_path: Path, grammar: str, headers: list[bytes]
) -> None:
    ref = _reference(tmp_path, headers)
    assert check_reference(str(ref), grammar) == []


def test_crlf_reference_is_clean(tmp_path: Path) -> None:
    # [S115]: a carriage return ending the line is allowed.
    ref = tmp_path / "crlf.fas"
    ref.write_bytes(
        STAMPA_OK + b"\r\n" + SEQ + b"\r\n"
        + b">refB Bacteria|Firmicutes\r\n" + SEQ + b"\r\n"
    )
    assert check_reference(str(ref), "stampa") == []


# ---------------------------------------------------------------------------
# errors ([S115]), always on a header other than the first
# ---------------------------------------------------------------------------

@pytest.mark.parametrize(
    ("grammar", "bad_header", "issue"),
    [
        ("stampa", b">refB", "format"),
        ("stampa", b">refB;tax=d:Bacteria,p:Firmicutes;", "format"),
        ("sintax", b">refB Bacteria|Firmicutes", "format"),
        ("sintax", b">refB;tax=;", "format"),
        ("sintax", b">refB;tax=d:Bacteria,Firmicutes;", "format"),
        ("stampa", b">refB Bacteria|Genus_\xe9", "invalid_utf8"),
        ("sintax", b">refB;tax=d:Bacteria,g:Genus_\xe9;", "invalid_utf8"),
        ("stampa", b">refB\tBacteria|Firmicutes", "control_char"),
        ("stampa", b">refB Bacteria|Firm\x01icutes", "control_char"),
        ("sintax", b">refB;tax=d:Bacteria,p:Firm\x1bicutes;", "control_char"),
        ("stampa", b">refA Bacteria|Firmicutes", "duplicate_accession"),
        ("sintax", b">refA;tax=d:Bacteria,p:Firmicutes;",
         "duplicate_accession"),
        ("stampa", b">refB  Bacteria|Firmicutes", "whitespace"),
        ("stampa", b">refB Bacteria|Firmicutes ", "whitespace"),
    ],
)
def test_error_on_a_later_header_is_reported(
    tmp_path: Path, grammar: str, bad_header: bytes, issue: str
) -> None:
    good = STAMPA_OK if grammar == "stampa" else SINTAX_OK
    ref = _reference(tmp_path, [good, bad_header])
    findings = check_reference(str(ref), grammar)
    assert issue in _issues(findings)
    # the finding points at the bad header: line 3 (header, seq, header)
    assert [f.line for f in findings if f.issue == issue] == [3]
    assert status(findings) == "error"


def test_duplicate_accession_is_reported_on_the_later_header(
    tmp_path: Path,
) -> None:
    ref = _reference(
        tmp_path,
        [STAMPA_OK, b">refB Bacteria|Firmicutes", b">refA Bacteria|Other"],
    )
    assert check_reference(str(ref), "stampa") == [
        Finding(5, "refA", "duplicate_accession")
    ]


def test_comma_in_accession_only_with_majority(tmp_path: Path) -> None:
    # [S115]: [S111] extended to the whole file, and only when the
    # majority step will look accessions up.
    hybrid = b">refB;tax=d:Bacteria,p:Firmicutes; Bacteria|Firmicutes"
    ref = _reference(tmp_path, [STAMPA_OK, hybrid])
    assert check_reference(str(ref), "stampa") == []
    assert check_reference(str(ref), "stampa", majority=True) == [
        Finding(3, "refB;tax=d:Bacteria,p:Firmicutes;", "comma_in_accession")
    ]


# ---------------------------------------------------------------------------
# warnings ([S116])
# ---------------------------------------------------------------------------

@pytest.mark.parametrize(
    ("grammar", "header", "issue"),
    [
        ("stampa", b">refB Bacteria||Gammaproteobacteria", "empty_rank"),
        ("stampa", b">refB Bacteria|Proteobacteria|", "empty_rank"),
        ("sintax", b">refB;tax=d:Bacteria,p:,c:Gamma;", "empty_rank"),
        ("sintax", b">refB;tax=d:Bacteria,,c:Gamma;", "empty_rank"),
        ("stampa", b">refB Bacteria|*|Gammaproteobacteria", "star_in_lineage"),
        ("stampa", b">refB Bacteria|Proteo*bacteria", "star_in_lineage"),
    ],
)
def test_warning_on_a_later_header_is_reported(
    tmp_path: Path, grammar: str, header: bytes, issue: str
) -> None:
    good = STAMPA_OK if grammar == "stampa" else SINTAX_OK
    ref = _reference(tmp_path, [good, header])
    findings = check_reference(str(ref), grammar)
    assert findings == [Finding(3, "refB", issue)]
    assert status(findings) == "warning"


# ---------------------------------------------------------------------------
# compressed references ([S114])
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("opener", [gzip.open, bz2.open])
def test_compressed_reference_is_checked_the_same_way(
    tmp_path: Path, opener
) -> None:
    ref = tmp_path / "reference.fas.compressed"
    with opener(ref, "wb") as handle:
        handle.write(STAMPA_OK + b"\n" + SEQ + b"\n"
                     + b">refB\tBacteria\n" + SEQ + b"\n")
    assert check_reference(str(ref), "stampa") == [
        Finding(3, "refB", "control_char")
    ]


# ---------------------------------------------------------------------------
# status and summary ([S117])
# ---------------------------------------------------------------------------

def test_status() -> None:
    assert status([]) == "clean"
    assert status([Finding(3, "refB", "empty_rank")]) == "warning"
    assert status([Finding(3, "refB", "empty_rank"),
                   Finding(5, "refC", "format")]) == "error"


def test_summary_names_parameter_counts_headers_and_report() -> None:
    findings = [Finding(2 * i + 1, f"ref{i}", "whitespace")
                for i in range(1, 8)]
    findings.append(Finding(17, "refZ", "empty_rank"))
    message = summarise(findings, "reference_dataset", "logs/x.tsv")
    assert "reference_dataset" in message
    assert "whitespace: 7" in message
    assert "empty_rank: 1" in message
    assert "logs/x.tsv" in message
    # the first five offending headers, with line and accession
    for i in range(1, 6):
        assert f"line {2 * i + 1}: ref{i}" in message
    assert "ref6" not in message


def test_summary_is_empty_for_a_clean_reference() -> None:
    assert summarise([], "reference_dataset", "logs/x.tsv") == ""


# ---------------------------------------------------------------------------
# command line ([S117])
# ---------------------------------------------------------------------------

def _run(tmp_path: Path, ref: Path, *extra: str) -> tuple[Path, Path]:
    report = tmp_path / "report.tsv"
    summary = tmp_path / "summary.txt"
    result = subprocess.run(
        [sys.executable, str(BIN), "--grammar", "stampa",
         "--parameter", "reference_dataset",
         "--report", str(report), "--summary", str(summary),
         *extra, str(ref)],
        capture_output=True,
    )
    # errors are reported, not raised: the report must be published
    # even when the run then aborts
    assert result.returncode == 0, result.stderr
    return report, summary


def test_cli_writes_report_and_summary_and_exits_zero_on_errors(
    tmp_path: Path,
) -> None:
    ref = _reference(tmp_path, [STAMPA_OK, b">refB Bacteria ",
                                b">refC Bacteria||Gamma"])
    report, summary = _run(tmp_path, ref)
    assert report.read_bytes().splitlines() == [
        b"line\taccession\tissue",
        b"3\trefB\twhitespace",
        b"5\trefC\tempty_rank",
    ]
    lines = summary.read_text().splitlines()
    assert lines[0] == "error"
    assert "reference_dataset" in "\n".join(lines[1:])


def test_cli_clean_reference_gives_header_only_report(tmp_path: Path) -> None:
    ref = _reference(tmp_path, [STAMPA_OK])
    report, summary = _run(tmp_path, ref)
    assert report.read_bytes().splitlines() == [b"line\taccession\tissue"]
    assert summary.read_text().splitlines() == ["clean"]


def test_cli_majority_flag_enables_comma_check(tmp_path: Path) -> None:
    ref = _reference(tmp_path, [b">a,b Bacteria"])
    report, _ = _run(tmp_path, ref, "--majority")
    assert b"1\ta,b\tcomma_in_accession" in report.read_bytes()


def test_cli_report_keeps_non_utf8_accession_bytes(tmp_path: Path) -> None:
    # [S112] spirit: the report shows the accession as it is in the file.
    ref = _reference(tmp_path, [b">ref\xe9 Bacteria"])
    report, _ = _run(tmp_path, ref)
    assert b"1\tref\xe9\tinvalid_utf8" in report.read_bytes()
