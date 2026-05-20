"""Unit tests for bin/stampa_merge.py."""

# COVERAGE: [S49]

from __future__ import annotations

import io
from pathlib import Path

import pytest

from stampa_merge import (
    emit_results,
    last_common_ancestor,
    parse_hits,
)


# ---------------------------------------------------------------------------
# Pure helpers
# ---------------------------------------------------------------------------

def test_last_common_ancestor_single_hit_passthrough() -> None:
    assert last_common_ancestor([["A", "B", "C"]]) == ["A", "B", "C"]


def test_last_common_ancestor_identical_lineages() -> None:
    # All hits agree → every level survives.
    taxonomies = [["A", "B", "C"], ["A", "B", "C"]]
    assert last_common_ancestor(taxonomies) == ["A", "B", "C"]


def test_last_common_ancestor_disagreeing_levels_become_star() -> None:
    # The phylum-level disagrees, so level 1 becomes "*"; level 0
    # still agrees → "A".
    taxonomies = [["A", "B", "C"], ["A", "X", "C"]]
    assert last_common_ancestor(taxonomies) == ["A", "*", "C"]


def test_last_common_ancestor_empty_returns_empty() -> None:
    assert last_common_ancestor([]) == []


# ---------------------------------------------------------------------------
# parse_hits — input shape
# ---------------------------------------------------------------------------

def _write(tmp_path: Path, name: str, content: str) -> Path:
    path = tmp_path / name
    path.write_text(content)
    return path


def test_parse_hits_strips_abundance_suffix(tmp_path: Path) -> None:
    hits = _write(tmp_path, "hits", "amp1_10\t99.5\trefA Kingdom|Phylum\n")
    rows = list(parse_hits(str(hits)))
    assert len(rows) == 1
    amplicon, abundance, identity, taxonomy, accession = rows[0]
    assert amplicon == "amp1"
    assert abundance == "10"
    assert identity == "99.5"
    assert taxonomy == ["Kingdom", "Phylum"]
    assert accession == "refA"


def test_parse_hits_handles_no_hit_marker(tmp_path: Path) -> None:
    hits = _write(tmp_path, "hits", "amp2_5\t0.0\t*\n")
    rows = list(parse_hits(str(hits)))
    assert rows == [("amp2", "5", "0.0", ["No_hit"], "No_hit")]


def test_parse_hits_handles_multi_underscore_amplicon(tmp_path: Path) -> None:
    # rsplit("_", 1) splits only on the final underscore so amplicon
    # IDs may themselves contain underscores (e.g. SHA1-derived
    # sample-prefixed IDs). The slurm-era stampa_merge.py used
    # plain split("_") and crashed on these.
    hits = _write(
        tmp_path, "hits",
        "sample_a_amp7_42\t99.0\trefX Bacteria|Firmicutes\n",
    )
    rows = list(parse_hits(str(hits)))
    assert rows == [
        ("sample_a_amp7", "42", "99.0", ["Bacteria", "Firmicutes"], "refX")
    ]


def test_parse_hits_malformed_hit_raises(tmp_path: Path) -> None:
    # Fail-loud on malformed input: a line that doesn't split into
    # exactly 3 tab-separated fields raises ValueError so upstream
    # silent data corruption is impossible to miss.
    hits = _write(tmp_path, "hits", "only\ttwo_columns\n")
    with pytest.raises(ValueError):
        list(parse_hits(str(hits)))


# ---------------------------------------------------------------------------
# End-to-end: emit_results
# ---------------------------------------------------------------------------

def test_emit_results_one_amplicon_one_hit(tmp_path: Path) -> None:
    hits = _write(
        tmp_path, "hits",
        "amp1_10\t99.5\trefA Kingdom|Phylum|Class\n",
    )
    buf = io.StringIO()
    emit_results(str(hits), buf)
    assert buf.getvalue() == (
        "amp1\t10\t99.5\tKingdom|Phylum|Class\trefA\n"
    )


def test_emit_results_groups_multi_hit_amplicon(tmp_path: Path) -> None:
    # Two top hits for amp1: they agree at level 0 and 2, disagree
    # at level 1 → LCA = "K|*|C".
    hits_text = (
        "amp1_10\t99.5\trefA K|P1|C\n"
        "amp1_10\t99.5\trefB K|P2|C\n"
    )
    hits = _write(tmp_path, "hits", hits_text)
    buf = io.StringIO()
    emit_results(str(hits), buf)
    assert buf.getvalue() == "amp1\t10\t99.5\tK|*|C\trefA,refB\n"


def test_emit_results_multiple_amplicons(tmp_path: Path) -> None:
    hits_text = (
        "amp1_10\t99.5\trefA K|P\n"
        "amp2_5\t88.0\t*\n"
    )
    hits = _write(tmp_path, "hits", hits_text)
    buf = io.StringIO()
    emit_results(str(hits), buf)
    assert buf.getvalue().splitlines() == [
        "amp1\t10\t99.5\tK|P\trefA",
        "amp2\t5\t88.0\tNo_hit\tNo_hit",
    ]


def test_emit_results_empty_input_emits_nothing(tmp_path: Path) -> None:
    hits = _write(tmp_path, "hits", "")
    buf = io.StringIO()
    emit_results(str(hits), buf)
    assert buf.getvalue() == ""
