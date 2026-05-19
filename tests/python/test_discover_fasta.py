"""Unit tests for ``bin/discover_fasta.py``.

Part B builds its fasta channel from the per-sample ``.fas`` files
produced by Part A (or, when run standalone, from the list of
folders passed via ``--fasta_folder``). The discoverer enforces two
rules:

* fasta files whose basename ends in ``_notmerged.fas`` are excluded
  — they belong to the shadow pipeline (``[S04]``) and are processed
  by a dedicated downstream path, not by the regular Part B
  pipeline;
* duplicate sample IDs (the basename minus the ``.fas`` extension)
  are an error and the workflow aborts.

See ``[S13]``, ``[S14]``, ``[S27]`` in SPECIFICATIONS.md.
"""

# COVERAGE: [S13], [S14], [S27]

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

from discover_fasta import (
    DuplicateSampleIDError,
    FastaSample,
    check_unique_sample_ids,
    discover,
)


def _make_fasta(folder: Path, name: str, body: str = ">x\nACGT\n") -> Path:
    p = folder / name
    p.write_text(body)
    return p


# ---------- discover() -----------------------------------------------------

def test_discover_lists_every_fasta_in_a_single_folder(tmp_path: Path) -> None:
    a = _make_fasta(tmp_path, "A.fas")
    b = _make_fasta(tmp_path, "B.fas")

    result = discover([tmp_path])

    assert result == [
        FastaSample(sample_id="A", fasta=a),
        FastaSample(sample_id="B", fasta=b),
    ]


def test_discover_excludes_notmerged_artefacts(tmp_path: Path) -> None:
    # [S27] — the shadow pipeline publishes <sampleId>_notmerged.fas
    # alongside the regular .fas. Part B's channel must skip them.
    a = _make_fasta(tmp_path, "A.fas")
    _make_fasta(tmp_path, "A_notmerged.fas")
    b = _make_fasta(tmp_path, "B.fas")

    result = discover([tmp_path])

    assert result == [
        FastaSample(sample_id="A", fasta=a),
        FastaSample(sample_id="B", fasta=b),
    ]


def test_discover_walks_every_folder_in_the_list(tmp_path: Path) -> None:
    # Same multi-folder semantics as discover_fastq.discover ([S10]).
    folder_a = tmp_path / "a"
    folder_b = tmp_path / "b"
    folder_a.mkdir()
    folder_b.mkdir()
    fa = _make_fasta(folder_a, "From_A.fas")
    fb = _make_fasta(folder_b, "from_b.fas")

    result = {s.sample_id: s for s in discover([folder_a, folder_b])}

    assert result == {
        "From_A": FastaSample(sample_id="From_A", fasta=fa),
        "from_b": FastaSample(sample_id="from_b", fasta=fb),
    }


def test_discover_ignores_non_fasta_files(tmp_path: Path) -> None:
    _make_fasta(tmp_path, "real.fas")
    (tmp_path / "real.qual").write_text("anything")
    (tmp_path / "real.stats").write_text("anything")
    (tmp_path / "notes.txt").write_text("ignore")

    result = discover([tmp_path])

    assert [s.fasta.name for s in result] == ["real.fas"]


def test_discover_keeps_empty_fasta_files(tmp_path: Path) -> None:
    # [S09] / [S27]: empty samples must travel through to the
    # occurrence table — they contribute a zero-filled sample column
    # (never a row). The fasta channel therefore keeps every .fas
    # regardless of size; downstream processes (build_distribution_file,
    # global_dereplication, etc.) are responsible for tolerating
    # empty inputs.
    real = _make_fasta(tmp_path, "real.fas")
    empty = tmp_path / "empty.fas"
    empty.write_text("")

    result = discover([tmp_path])

    assert result == [
        FastaSample(sample_id="empty", fasta=empty),
        FastaSample(sample_id="real", fasta=real),
    ]


# ---------- uniqueness ([S13]/[S14]) --------------------------------------

def test_check_unique_sample_ids_passes_on_distinct_ids() -> None:
    samples = [
        FastaSample(sample_id="A", fasta=Path("/run1/A.fas")),
        FastaSample(sample_id="B", fasta=Path("/run1/B.fas")),
    ]
    check_unique_sample_ids(samples)


def test_check_unique_sample_ids_raises_on_duplicates() -> None:
    samples = [
        FastaSample(sample_id="A", fasta=Path("/run1/A.fas")),
        FastaSample(sample_id="A", fasta=Path("/run2/A.fas")),
        FastaSample(sample_id="B", fasta=Path("/run1/B.fas")),
    ]
    with pytest.raises(DuplicateSampleIDError) as excinfo:
        check_unique_sample_ids(samples)
    msg = str(excinfo.value)
    assert "A" in msg
    assert "/run1/A.fas" in msg
    assert "/run2/A.fas" in msg
    # the non-duplicated sample's path is not surfaced
    assert "/run1/B.fas" not in msg


def test_cli_emits_tsv_for_each_sample(tmp_path: Path) -> None:
    a = _make_fasta(tmp_path, "A.fas")
    b = _make_fasta(tmp_path, "B.fas")
    # _notmerged is dropped
    _make_fasta(tmp_path, "A_notmerged.fas")

    script = (
        Path(__file__).resolve().parents[2] / "bin" / "discover_fasta.py"
    )
    proc = subprocess.run(
        [sys.executable, str(script), str(tmp_path)],
        capture_output=True,
        text=True,
    )
    assert proc.returncode == 0, proc.stderr
    lines = sorted(proc.stdout.strip().splitlines())
    assert lines == [
        f"A\t{a}",
        f"B\t{b}",
    ]


def test_cli_exits_non_zero_on_duplicate_sample_ids(tmp_path: Path) -> None:
    a = tmp_path / "run1"
    b = tmp_path / "run2"
    a.mkdir()
    b.mkdir()
    _make_fasta(a, "A.fas")
    _make_fasta(b, "A.fas")

    script = (
        Path(__file__).resolve().parents[2] / "bin" / "discover_fasta.py"
    )
    proc = subprocess.run(
        [sys.executable, str(script), str(a), str(b)],
        capture_output=True,
        text=True,
    )
    assert proc.returncode != 0, (
        f"expected non-zero exit; stdout={proc.stdout!r}"
    )
    assert "duplicate" in proc.stderr.lower(), proc.stderr
    assert str(a / "A.fas") in proc.stderr
    assert str(b / "A.fas") in proc.stderr
