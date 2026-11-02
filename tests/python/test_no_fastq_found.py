"""CLI tests for ``bin/discover_fastq.py`` when the scan finds no
fastq file at all ([S126]).

Discovery is not recursive ([S10]/[S21]): a ``--fastq_folder`` that
points at the parent of the real input folders used to yield an empty
samples list, and Part A then finished "successfully" having processed
nothing. The helper must instead exit non-zero with a message that
names the scanned folder(s), says discovery is not recursive, and
points at any subfolder that does hold fastq files.
"""

# COVERAGE: [S126]

import subprocess
import sys
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[2] / "bin" / "discover_fastq.py"


def _run(*folders: Path) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, str(SCRIPT), *[str(f) for f in folders]],
        capture_output=True,
        text=True,
    )


def _touch(path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("")
    return path


def test_empty_folder_exits_non_zero(tmp_path: Path) -> None:
    proc = _run(tmp_path)

    assert proc.returncode != 0, f"expected failure; stdout={proc.stdout!r}"
    assert proc.stdout == ""
    assert str(tmp_path) in proc.stderr, proc.stderr
    assert "not recursive" in proc.stderr.lower(), proc.stderr


def test_folder_with_only_non_fastq_files_counts_as_empty(
    tmp_path: Path,
) -> None:
    for name in ["S.fas", "S.qual", "S.stats", "S.log"]:
        _touch(tmp_path / name)

    proc = _run(tmp_path)

    assert proc.returncode != 0, f"expected failure; stdout={proc.stdout!r}"
    assert proc.stdout == ""
    assert "not recursive" in proc.stderr.lower(), proc.stderr


def test_subfolders_holding_fastq_files_are_listed(tmp_path: Path) -> None:
    # the layout that triggered [S126]: inputs one level down
    for sub in ["PNA", "nonPNA"]:
        _touch(tmp_path / sub / f"S_{sub}.R1.fq.gz")
        _touch(tmp_path / sub / f"S_{sub}.R2.fq.gz")
        _touch(tmp_path / sub / f"S_{sub}.fas")
    _touch(tmp_path / "unrelated" / "notes.txt")

    proc = _run(tmp_path)

    assert proc.returncode != 0, f"expected failure; stdout={proc.stdout!r}"
    assert proc.stdout == ""
    assert str(tmp_path / "PNA") in proc.stderr, proc.stderr
    assert str(tmp_path / "nonPNA") in proc.stderr, proc.stderr
    # a subfolder without fastq files is not suggested
    assert str(tmp_path / "unrelated") not in proc.stderr, proc.stderr


def test_every_scanned_folder_is_named(tmp_path: Path) -> None:
    first = tmp_path / "run1"
    second = tmp_path / "run2"
    first.mkdir()
    second.mkdir()

    proc = _run(first, second)

    assert proc.returncode != 0, f"expected failure; stdout={proc.stdout!r}"
    assert str(first) in proc.stderr, proc.stderr
    assert str(second) in proc.stderr, proc.stderr


def test_partially_empty_scan_warns_about_the_empty_folders(
    tmp_path: Path,
) -> None:
    # [S126] fails only when the whole scan is empty; an empty folder
    # next to a populated one is reported by a warning, not an error.
    empty = tmp_path / "empty"
    empty.mkdir()
    fastq = _touch(tmp_path / "full" / "loose_sample.fastq.gz")

    proc = _run(empty, fastq.parent)

    assert proc.returncode == 0, proc.stderr
    assert proc.stdout.startswith("loose_sample\t"), proc.stdout
    warnings = [
        line for line in proc.stderr.splitlines()
        if line.startswith("warning:")
    ]
    assert len(warnings) == 1, proc.stderr
    assert str(empty) in warnings[0], proc.stderr
    assert str(fastq.parent) not in proc.stderr, proc.stderr


def test_populated_scan_prints_no_warning(tmp_path: Path) -> None:
    _touch(tmp_path / "loose_sample.fastq.gz")

    proc = _run(tmp_path)

    assert proc.returncode == 0, proc.stderr
    assert proc.stderr == "", proc.stderr
