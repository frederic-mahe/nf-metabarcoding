"""Characterization tests for ``bin/merge_substring_otus.py``.

Byte-exact golden-file gate that pins the legacy script's current
behaviour before any refactor. The main fixture exercises:

* a pass-through OTU (neither master nor pupil) → emitted as-is
* a master with a single pupil → samples summed onto the master,
  ``spread`` recomputed from the merged columns, ``total`` summed,
  ``cloud`` incremented by ``pupil_cloud + 1`` (the "+1 per merged
  pupil" quirk)
* a master with two pupils → incremental accumulation holds
* output order: pass-through rows in input order, then master rows
  in insertion order (the bash sort step downstream restores the
  numeric OTU ordering)

A separate test exercises the overlap branch: when an OTU is both
master and pupil in the match list, the script must exit non-zero
and print the documented WARNING line on stderr.

The match file format mirrors what vsearch ``--uc`` emits as ``^H``
lines: 10 tab-separated fields ending in ``pupil\\tmaster``; the
script reads ``line.strip().split("\\t")[-2:]``.
"""

# COVERAGE: [S39]

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


SCRIPT = (
    Path(__file__).resolve().parents[2]
    / "bin"
    / "merge_substring_otus.py"
)
FIXTURE = (
    Path(__file__).resolve().parents[1]
    / "data"
    / "merge_substring_otus"
)


def test_golden_output(tmp_path: Path) -> None:
    out = tmp_path / "merged_table"
    subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "-t", str(FIXTURE / "otu_table"),
            "-m", str(FIXTURE / "matches"),
            "-o", str(out),
        ],
        check=True,
    )
    expected = (FIXTURE / "expected_output").read_text()
    assert out.read_text() == expected, (
        "byte-exact mismatch.\n"
        f"Got:\n{out.read_text()!r}\n"
        f"Expected:\n{expected!r}"
    )


def test_master_pupil_overlap_aborts(tmp_path: Path) -> None:
    out = tmp_path / "out"
    proc = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "-t", str(FIXTURE / "otu_table"),
            "-m", str(FIXTURE / "matches_overlap"),
            "-o", str(out),
        ],
        capture_output=True,
        text=True,
    )
    assert proc.returncode != 0, (
        f"expected non-zero exit; stdout={proc.stdout!r}"
    )
    assert "WARNING" in proc.stderr
    assert "common to hit and query columns" in proc.stderr
