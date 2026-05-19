"""Characterization tests for ``tmp/stampa/stampa_merge.py``.

These tests lock the legacy script's behaviour as we found it
(slurm-array shard loop + last-common-ancestor + empty-file
abort). They are the safety net for the refactor toward
``bin/stampa_merge.py`` described in [S49] — see TODO.md / chat
history (Plan A option c, 2026-05-19).

Two strategies are used:

* import-based tests for the pure ``last_common_ancestor`` helper,
  importing the legacy module via a shim that suppresses the
  module-level ``sys.exit(0)`` (which would otherwise terminate
  the pytest process when the module is imported).
* subprocess-based tests for ``main()``, since it does
  ``sys.argv[1]`` + ``os.chdir`` + ``sys.exit`` and writes files
  into the working directory — a fully isolated process is the
  honest interface.
"""

# COVERAGE: [S49]

from __future__ import annotations

import importlib.util
import subprocess
import sys
import types
from pathlib import Path
from unittest.mock import patch

import pytest


# Path to the legacy script under test. ``tests/python/...`` ->
# repo-root walk: parent (python) -> parent (tests) -> parent (root).
LEGACY_PATH: Path = (
    Path(__file__).resolve().parent.parent.parent
    / "tmp" / "stampa" / "stampa_merge.py"
)


def _import_legacy() -> types.ModuleType:
    """Import the legacy module while suppressing its top-level sys.exit(0).

    The legacy script ends with a bare ``sys.exit(0)`` outside the
    ``if __name__ == '__main__'`` guard. Patching ``sys.exit`` to a
    no-op for the duration of ``exec_module`` lets us load the module
    for unit testing without terminating the pytest process.
    """
    spec = importlib.util.spec_from_file_location(
        "legacy_stampa_merge", LEGACY_PATH
    )
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    with patch.object(sys, "exit", lambda *args, **kwargs: None):
        spec.loader.exec_module(mod)
    return mod


legacy = _import_legacy()


# ---------------------------------------------------------------------------
# last_common_ancestor — pure function
# ---------------------------------------------------------------------------

def test_lca_single_taxonomy_passthrough() -> None:
    """A single top-hit's taxonomy is returned verbatim."""
    assert legacy.last_common_ancestor([["A", "B", "C"]]) == ["A", "B", "C"]


def test_lca_all_taxonomies_agree() -> None:
    """Every top hit shares the same taxonomy → all levels survive."""
    taxonomies = [["A", "B", "C"], ["A", "B", "C"]]
    assert legacy.last_common_ancestor(taxonomies) == ["A", "B", "C"]


def test_lca_disagreement_rewrites_level_to_star() -> None:
    """A divergent level becomes ``*``; matching levels stay."""
    taxonomies = [["A", "B", "C"], ["A", "X", "C"]]
    assert legacy.last_common_ancestor(taxonomies) == ["A", "*", "C"]


def test_lca_multiple_disagreements() -> None:
    """Each divergent level independently becomes ``*``."""
    taxonomies = [
        ["Bacteria", "Firmicutes", "Bacilli"],
        ["Bacteria", "Proteobacteria", "Bacilli"],
        ["Bacteria", "Firmicutes", "Clostridia"],
    ]
    assert legacy.last_common_ancestor(taxonomies) == ["Bacteria", "*", "*"]


def test_lca_empty_input_raises_index_error() -> None:
    """Legacy quirk: an empty input list raises ``IndexError``.

    The script's ``else`` branch indexes ``taxonomies[0]`` and the
    except clause re-raises after writing the input to stderr. Pinned
    here so the refactor knows this code path exists — the new
    ``bin/stampa_merge.py`` returns ``[]`` instead, a deliberate
    behaviour change that is *not* a regression.
    """
    with pytest.raises(IndexError):
        legacy.last_common_ancestor([])


# ---------------------------------------------------------------------------
# main — subprocess end-to-end
# ---------------------------------------------------------------------------

def _run_legacy(directory: Path) -> subprocess.CompletedProcess[str]:
    """Run the legacy script against *directory* in a child process."""
    return subprocess.run(
        [sys.executable, str(LEGACY_PATH), str(directory)],
        capture_output=True,
        text=True,
    )


def test_main_writes_results_per_shard(tmp_path: Path) -> None:
    """Each ``hits.NNNNN`` produces a sibling ``results.NNNNN``."""
    (tmp_path / "hits.00000").write_text(
        "amp1_10\t99.0\trefA Bacteria|Firmicutes\n"
        "amp1_10\t99.0\trefB Bacteria|Firmicutes\n"
        "amp2_5\t95.0\trefC Bacteria|Proteobacteria\n"
    )
    (tmp_path / "hits.00001").write_text(
        "amp3_2\t100.0\trefD Eukaryota|Fungi\n"
    )

    result = _run_legacy(tmp_path)

    assert result.returncode == 0, result.stderr
    assert (tmp_path / "results.00000").read_text() == (
        "amp1\t10\t99.0\tBacteria|Firmicutes\trefA,refB\n"
        "amp2\t5\t95.0\tBacteria|Proteobacteria\trefC\n"
    )
    assert (tmp_path / "results.00001").read_text() == (
        "amp3\t2\t100.0\tEukaryota|Fungi\trefD\n"
    )


def test_main_handles_no_hit_marker(tmp_path: Path) -> None:
    """A ``*`` hit field rewrites accession + taxonomy to ``No_hit``."""
    (tmp_path / "hits.00000").write_text("amp1_3\t0.0\t*\n")

    result = _run_legacy(tmp_path)

    assert result.returncode == 0, result.stderr
    assert (tmp_path / "results.00000").read_text() == (
        "amp1\t3\t0.0\tNo_hit\tNo_hit\n"
    )


def test_main_lca_collapses_divergent_levels(tmp_path: Path) -> None:
    """Two hits with divergent phyla → ``*`` at that level in the output."""
    (tmp_path / "hits.00000").write_text(
        "amp1_4\t99.0\trefA Bacteria|Firmicutes|Bacilli\n"
        "amp1_4\t99.0\trefB Bacteria|Proteobacteria|Bacilli\n"
    )

    result = _run_legacy(tmp_path)

    assert result.returncode == 0, result.stderr
    assert (tmp_path / "results.00000").read_text() == (
        "amp1\t4\t99.0\tBacteria|*|Bacilli\trefA,refB\n"
    )


def test_main_empty_shard_exits_non_zero(tmp_path: Path) -> None:
    """A zero-byte ``hits.*`` aborts the run with stderr listing it."""
    (tmp_path / "hits.00000").write_text("amp1_1\t100.0\t*\n")
    (tmp_path / "hits.00001").write_text("")  # legitimate failure mode

    result = _run_legacy(tmp_path)

    assert result.returncode == 1
    assert "empty hits file" in result.stderr
    assert "hits.00001" in result.stderr


def test_main_missing_directory_exits_non_zero(tmp_path: Path) -> None:
    """Pointing at a non-existent directory exits with a clear error."""
    target = tmp_path / "does-not-exist"

    result = _run_legacy(target)

    assert result.returncode != 0
    # Legacy uses ``sys.exit("ERROR: directory ... not found!")``,
    # which writes the message to stderr.
    assert "not found" in result.stderr


def test_main_shard_walk_is_sorted(tmp_path: Path) -> None:
    """``hits.*`` files are processed in alphabetical order.

    Pinning this so the refactor doesn't accidentally pick a
    different ordering (e.g. inode order from ``os.scandir``).
    Amplicon IDs deliberately use exactly one underscore (the
    ``;size=`` → ``_`` rewrite shape from ``stampa.sh``) since the
    legacy split-on-underscore crashes otherwise — see
    ``test_main_multi_underscore_amplicon_is_a_known_limitation``.
    """
    (tmp_path / "hits.00002").write_text("ampC_1\t90.0\t*\n")
    (tmp_path / "hits.00000").write_text("ampA_1\t90.0\t*\n")
    (tmp_path / "hits.00001").write_text("ampB_1\t90.0\t*\n")

    result = _run_legacy(tmp_path)

    assert result.returncode == 0, result.stderr
    # Per-shard content is independent; we just check each shard's
    # output begins with its expected amplicon prefix.
    for index, expected_prefix in enumerate(["ampA", "ampB", "ampC"]):
        content = (tmp_path / f"results.0000{index}").read_text()
        assert content.startswith(expected_prefix), (
            f"results.0000{index} unexpected content: {content!r}"
        )


def test_main_malformed_hit_aborts(tmp_path: Path) -> None:
    """A hit line that doesn't split into 3 tab-separated fields aborts.

    Legacy does ``amplicon, identity, hit = line.strip().split('\\t')``
    with no except, so any short line raises ``ValueError`` and the
    process exits non-zero. Pinning here so the refactor preserves
    this "fail loud on malformed input" property.
    """
    (tmp_path / "hits.00000").write_text("only\ttwo_columns\n")

    result = _run_legacy(tmp_path)

    assert result.returncode != 0


def test_main_multi_underscore_amplicon_is_a_known_limitation(
    tmp_path: Path,
) -> None:
    """Legacy gap: amplicon IDs with multiple underscores crash on split.

    ``amp_foo_10`` splits into 3 fields → ``ValueError`` on the
    ``amplicon, abundance = amplicon.split('_')`` unpacking. Pinned
    so the refactor knows to switch to ``rsplit('_', 1)`` (which the
    new ``bin/stampa_merge.py`` already does, see [S49]) and flip
    this test to assert success.
    """
    (tmp_path / "hits.00000").write_text("amp_foo_10\t90.0\t*\n")

    result = _run_legacy(tmp_path)

    assert result.returncode != 0
