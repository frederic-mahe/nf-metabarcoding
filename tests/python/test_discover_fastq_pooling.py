"""CLI tests for the Nanopore-oriented options of ``bin/discover_fastq.py``.

* ``--pool-by-folder`` ([S128]): each first-level subfolder of a listed
  folder becomes one single-end sample named after it; the helper emits
  one ``sample<TAB>file<TAB>`` row per file.
* the numbered-series warning ([S130]): without pooling, a folder
  holding ``<prefix>_0``, ``_1``, ``_2`` ... triggers a warning.
* ``--recursive`` ([S131]): the scan descends into every subfolder.
"""

import subprocess
import sys
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[2] / "bin" / "discover_fastq.py"


def _run(*args: str) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, str(SCRIPT), *args],
        capture_output=True,
        text=True,
    )


def _touch(path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("")
    return path


def _rows(stdout: str) -> list[list[str]]:
    return [line.split("\t") for line in stdout.splitlines()]


def _nanopore_layout(root: Path) -> Path:
    """Build a MinKNOW-like fastq_pass/ folder under ``root``."""
    top = root / "fastq_pass"
    for n in range(3):
        _touch(top / "barcode01" / f"run_barcode01_{n}.fastq.gz")
    _touch(top / "barcode02" / "run_barcode02_0.fastq")
    _touch(top / "unclassified" / "run_unclassified_0.fastq.bz2")
    _touch(top / "other" / "notes.txt")
    _touch(top / "loose.fastq.gz")
    return top


# ---------- [S128] --pool-by-folder ----------------------------------------

# COVERAGE: [S128]
def test_pool_by_folder_one_sample_per_subfolder(tmp_path: Path) -> None:
    top = _nanopore_layout(tmp_path)

    proc = _run("--pool-by-folder", str(top))

    assert proc.returncode == 0, proc.stderr
    rows = _rows(proc.stdout)
    assert [r[0] for r in rows] == [
        "barcode01", "barcode01", "barcode01", "barcode02", "unclassified",
    ], rows
    # one file per row, sorted by path, no R2
    assert [Path(r[1]).name for r in rows] == [
        "run_barcode01_0.fastq.gz",
        "run_barcode01_1.fastq.gz",
        "run_barcode01_2.fastq.gz",
        "run_barcode02_0.fastq",
        "run_unclassified_0.fastq.bz2",
    ], rows
    assert all(len(r) == 3 and r[2] == "" for r in rows), rows


# COVERAGE: [S128]
def test_pool_by_folder_warns_about_loose_files_and_empty_subfolders(
    tmp_path: Path,
) -> None:
    top = _nanopore_layout(tmp_path)

    proc = _run("--pool-by-folder", str(top))

    warnings = [
        line for line in proc.stderr.splitlines()
        if line.startswith("warning: ")
    ]
    assert len(warnings) == 2, proc.stderr
    assert any(
        str(top / "other") in w for w in warnings
    ), proc.stderr
    assert any(
        str(top) in w and "1 " in w and "ignored" in w for w in warnings
    ), proc.stderr
    # no chunk-series warning under pooling ([S130])
    assert "numbered file series" not in proc.stderr, proc.stderr


# COVERAGE: [S128]
def test_pool_by_folder_without_fastq_subfolders_fails(
    tmp_path: Path,
) -> None:
    _touch(tmp_path / "loose.fastq.gz")
    _touch(tmp_path / "other" / "notes.txt")

    proc = _run("--pool-by-folder", str(tmp_path))

    assert proc.returncode != 0, proc.stdout
    assert proc.stdout == ""
    assert str(tmp_path) in proc.stderr, proc.stderr


# COVERAGE: [S128]
def test_pool_by_folder_duplicate_folder_names_fail(tmp_path: Path) -> None:
    for run in ("run1", "run2"):
        _touch(tmp_path / run / "barcode01" / "x_0.fastq.gz")

    proc = _run(
        "--pool-by-folder", str(tmp_path / "run1"), str(tmp_path / "run2"),
    )

    assert proc.returncode != 0, proc.stdout
    assert "barcode01" in proc.stderr, proc.stderr
    assert "duplicate" in proc.stderr.lower(), proc.stderr


# COVERAGE: [S128]
def test_pool_by_folder_checks_sample_id_charset(tmp_path: Path) -> None:
    _touch(tmp_path / "bad name" / "x_0.fastq.gz")

    proc = _run("--pool-by-folder", str(tmp_path))

    assert proc.returncode != 0, proc.stdout
    assert "bad name" in proc.stderr, proc.stderr


# ---------- [S130] numbered-series warning ---------------------------------

# COVERAGE: [S130]
def test_numbered_series_warns(tmp_path: Path) -> None:
    for n in range(3):
        _touch(tmp_path / f"x_{n}.fastq.gz")

    proc = _run(str(tmp_path))

    assert proc.returncode == 0, proc.stderr
    warnings = [
        line for line in proc.stderr.splitlines()
        if line.startswith("warning: ") and "--pool_fastq_by_folder" in line
    ]
    assert len(warnings) == 1, proc.stderr
    assert str(tmp_path) in warnings[0], warnings
    assert "x" in warnings[0], warnings


# COVERAGE: [S130]
def test_plain_pair_does_not_warn(tmp_path: Path) -> None:
    _touch(tmp_path / "s_1.fastq.gz")
    _touch(tmp_path / "s_2.fastq.gz")

    proc = _run(str(tmp_path))

    assert proc.returncode == 0, proc.stderr
    assert "warning" not in proc.stderr.lower(), proc.stderr


# COVERAGE: [S130]
def test_series_without_zero_does_not_warn(tmp_path: Path) -> None:
    for n in (1, 2, 3):
        _touch(tmp_path / f"x_{n}.fastq.gz")

    proc = _run(str(tmp_path))

    assert "--pool_fastq_by_folder" not in proc.stderr, proc.stderr


# ---------- [S131] --recursive ---------------------------------------------

# COVERAGE: [S131]
def test_recursive_finds_pair_two_levels_down(tmp_path: Path) -> None:
    deep = tmp_path / "a" / "b"
    _touch(deep / "S_R1.fastq.gz")
    _touch(deep / "S_R2.fastq.gz")

    proc = _run("--recursive", str(tmp_path))

    assert proc.returncode == 0, proc.stderr
    rows = _rows(proc.stdout)
    assert len(rows) == 1, rows
    sample, r1, r2 = rows[0]
    assert sample == "S"
    assert Path(r1) == deep / "S_R1.fastq.gz"
    assert Path(r2) == deep / "S_R2.fastq.gz"
    # the listed folder's tree holds fastq files: no empty-folder warning
    assert "warning" not in proc.stderr.lower(), proc.stderr


# COVERAGE: [S131]
def test_without_recursive_the_same_tree_fails(tmp_path: Path) -> None:
    deep = tmp_path / "a" / "b"
    _touch(deep / "S_R1.fastq.gz")
    _touch(deep / "S_R2.fastq.gz")

    proc = _run(str(tmp_path))

    assert proc.returncode != 0, proc.stdout
    assert "not recursive" in proc.stderr.lower(), proc.stderr


# COVERAGE: [S131]
def test_recursive_empty_tree_message(tmp_path: Path) -> None:
    _touch(tmp_path / "a" / "notes.txt")

    proc = _run("--recursive", str(tmp_path))

    assert proc.returncode != 0, proc.stdout
    assert "recursive" in proc.stderr.lower(), proc.stderr
    assert "not recursive" not in proc.stderr.lower(), proc.stderr


# COVERAGE: [S131]
def test_recursive_pairs_only_within_a_folder(tmp_path: Path) -> None:
    _touch(tmp_path / "a" / "S_R1.fastq.gz")
    _touch(tmp_path / "b" / "S_R2.fastq.gz")

    proc = _run("--recursive", str(tmp_path))

    assert proc.returncode == 0, proc.stderr
    rows = _rows(proc.stdout)
    # R1 has no partner in its own folder → single-end; the lone R2 is
    # single-end too; their derived IDs differ, so no collision.
    assert all(r[2] == "" for r in rows), rows
    assert len(rows) == 2, rows


# COVERAGE: [S131]
def test_recursive_collision_across_folders_fails(tmp_path: Path) -> None:
    for sub in ("a", "b"):
        _touch(tmp_path / sub / "S_R1.fastq.gz")
        _touch(tmp_path / sub / "S_R2.fastq.gz")

    proc = _run("--recursive", str(tmp_path))

    assert proc.returncode != 0, proc.stdout
    assert "duplicate" in proc.stderr.lower(), proc.stderr
