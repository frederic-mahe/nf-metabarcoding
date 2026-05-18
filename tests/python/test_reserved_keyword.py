"""Unit tests for the `notmerged` reserved-suffix rule in
``bin/discover_fastq.py`` ([S23]).

A user-supplied sample ID whose basename resolves to one ending in
the literal string ``notmerged`` must be rejected at discovery time,
before any Part A process runs. This keeps the shadow-pipeline
artefacts (`<sampleId>_notmerged.*`) collision-free with user
artefacts.
"""

# COVERAGE: [S23]

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

from discover_fastq import (
    check_reserved_suffix,
    ReservedSuffixError,
    Sample,
)


# ---------- pure-function check -------------------------------------------

def test_check_reserved_suffix_passes_on_clean_ids() -> None:
    samples = [
        Sample(
            sample_id="A",
            r1=Path("A_1.fastq.gz"),
            r2=Path("A_2.fastq.gz"),
        ),
        Sample(
            sample_id="B_run17",
            r1=Path("B_run17_1.fastq.gz"),
            r2=None,
        ),
    ]
    # no exception
    check_reserved_suffix(samples)


@pytest.mark.parametrize(
    "bad_id",
    [
        "X_notmerged",
        "notmerged",
        "foo_bar_notmerged",
        "SampleA_NOTMERGED".lower(),
    ],
)
def test_check_reserved_suffix_rejects_notmerged_suffix(bad_id: str) -> None:
    samples = [
        Sample(sample_id=bad_id, r1=Path(f"{bad_id}_1.fastq.gz"), r2=None),
    ]
    with pytest.raises(ReservedSuffixError) as excinfo:
        check_reserved_suffix(samples)
    assert "notmerged" in str(excinfo.value)
    assert bad_id in str(excinfo.value)


def test_check_reserved_suffix_reports_every_offender() -> None:
    samples = [
        Sample(
            sample_id="A_notmerged",
            r1=Path("A_notmerged_1.fastq.gz"),
            r2=None,
        ),
        Sample(sample_id="B", r1=Path("B_1.fastq.gz"), r2=None),
        Sample(
            sample_id="C_notmerged",
            r1=Path("C_notmerged_1.fastq.gz"),
            r2=None,
        ),
    ]
    with pytest.raises(ReservedSuffixError) as excinfo:
        check_reserved_suffix(samples)
    msg = str(excinfo.value)
    assert "A_notmerged" in msg
    assert "C_notmerged" in msg


# ---------- CLI behaviour --------------------------------------------------

def test_cli_exits_non_zero_when_reserved_suffix_present(
    tmp_path: Path,
) -> None:
    # Build a tiny folder with a paired-end pair whose R1 resolves to
    # `X_notmerged` via canonical pattern row 7 (`[._]1.<ext>`).
    r1 = tmp_path / "X_notmerged_1.fastq.gz"
    r2 = tmp_path / "X_notmerged_2.fastq.gz"
    r1.write_bytes(b"")
    r2.write_bytes(b"")

    script = (
        Path(__file__).resolve().parents[2] / "bin" / "discover_fastq.py"
    )
    proc = subprocess.run(
        [sys.executable, str(script), str(tmp_path)],
        capture_output=True,
        text=True,
    )
    assert proc.returncode != 0, (
        f"expected non-zero exit; stdout={proc.stdout!r}"
    )
    assert "notmerged" in proc.stderr.lower(), proc.stderr
