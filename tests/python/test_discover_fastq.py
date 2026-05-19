"""Unit tests for ``bin/discover_fastq.py``.

The pattern table is the single source of truth for both R2 name
derivation and sample-ID extraction. Cases here are ported from the
reference shell script
``/tmp/fred_01_merge_trim_and_dereplicate.sh`` (functions
``get_reverse_name`` and ``get_sample_name``), so this file is also
the spec-by-example for the canonical patterns documented in
``SPECIFICATIONS.md``.
"""

# COVERAGE: [S10], [S11], [S12]

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

from discover_fastq import (
    DuplicateSampleIDError,
    Sample,
    check_unique_sample_ids,
    derive_r2_name,
    derive_sample_id,
    discover,
)


# ---------- derive_r2_name -------------------------------------------------

# Every extension form must round-trip through pattern row 1.
@pytest.mark.parametrize(
    "r1, expected_r2",
    [
        ("A_L001_R1_001.fastq",         "A_L001_R2_001.fastq"),
        ("A_L001_R1_001.fastq.gz",      "A_L001_R2_001.fastq.gz"),
        ("A_L001_R1_001.fastq.bz2",     "A_L001_R2_001.fastq.bz2"),
        ("A_L001_R1_001.fq",            "A_L001_R2_001.fq"),
        ("A_L001_R1_001.fq.gz",         "A_L001_R2_001.fq.gz"),
        ("A_L001_R1_001.fq.bz2",        "A_L001_R2_001.fq.bz2"),
    ],
)
def test_derive_r2_name_accepts_every_extension(
    r1: str, expected_r2: str,
) -> None:
    assert derive_r2_name(r1) == expected_r2


@pytest.mark.parametrize(
    "r1, expected_r2",
    [
        # row 1: MiSeq w/ trailing 00x
        ("A_L001_R1_001.fastq.gz",      "A_L001_R2_001.fastq.gz"),
        ("A_L001_R1_002.fastq.gz",      "A_L001_R2_002.fastq.gz"),
        ("A_L009_R1_001.fastq.gz",      "A_L009_R2_001.fastq.gz"),
        # row 2: MiSeq w/ middle junk
        ("A_L001_junk_R1.fastq.gz",     "A_L001_junk_R2.fastq.gz"),
        ("A_L009_junk_R1.fastq.gz",     "A_L009_junk_R2.fastq.gz"),
        # row 3: plain MiSeq
        ("A_L001_R1.fastq.gz",          "A_L001_R2.fastq.gz"),
        ("A_L009_R1.fastq.gz",          "A_L009_R2.fastq.gz"),
        # row 4: numeric-lane w/ tail
        ("A_1_1_junk.fastq.gz",         "A_1_2_junk.fastq.gz"),
        ("A_9_1_junk.fastq.gz",         "A_9_2_junk.fastq.gz"),
        ("A.1_1_junk.fastq.gz",         "A.1_2_junk.fastq.gz"),
        # row 5: numeric-lane
        ("A_1_1.fastq.gz",              "A_1_2.fastq.gz"),
        ("A_9_1.fastq.gz",              "A_9_2.fastq.gz"),
        ("A.1_1.fastq.gz",              "A.1_2.fastq.gz"),
        # row 6: generic R1
        ("A_R1.fastq.gz",               "A_R2.fastq.gz"),
        ("A.R1.fastq.gz",               "A.R2.fastq.gz"),
    ],
)
def test_derive_r2_name_covers_every_canonical_row(
    r1: str, expected_r2: str,
) -> None:
    assert derive_r2_name(r1) == expected_r2


def test_derive_r2_name_returns_none_when_no_pattern_matches() -> None:
    # Ported from the bash script: `A_R3.fastq.gz` is not paired-end.
    assert derive_r2_name("A_R3.fastq.gz") is None


# ---------- derive_sample_id -----------------------------------------------

@pytest.mark.parametrize(
    "r1, expected_sample",
    [
        # row 1
        ("A_L001_R1_001.fastq.gz",       "A"),
        ("A_L009_R1_001.fastq.gz",       "A"),
        ("A_L001_R1_002.fastq.gz",       "A"),
        # row 2
        ("A_L001_junk_R1.fastq.gz",      "A"),
        ("A_L009_junk_R1.fastq.gz",      "A"),
        # row 3
        ("A_L001_R1.fastq.gz",           "A"),
        ("A_L009_R1.fastq.gz",           "A"),
        # row 4
        ("A_1_1_junk.fastq.gz",          "A"),
        ("A_9_1_junk.fastq.gz",          "A"),
        # row 5
        ("A_1_1.fastq.gz",               "A"),
        ("A_9_1.fastq.gz",               "A"),
        # row 6
        ("A_R1.fastq.gz",                "A"),
        ("A.R1.fastq.gz",                "A"),
        # row 7
        ("A_1.fastq.gz",                 "A"),
        ("A.1.fastq.gz",                 "A"),
    ],
)
def test_derive_sample_id_for_each_canonical_row(
    r1: str, expected_sample: str,
) -> None:
    assert derive_sample_id(r1) == expected_sample


@pytest.mark.parametrize(
    "r1, expected_sample",
    [
        # The "R1" or "1" token appears in the *sample* prefix as well —
        # the longest-suffix anchoring must pick the trailing token.
        ("R1_L001_R1_001.fastq.gz",      "R1"),
        ("A_R1_L001_R1_001.fastq.gz",    "A_R1"),
    ],
)
def test_derive_sample_id_anchors_on_trailing_pattern(
    r1: str, expected_sample: str,
) -> None:
    assert derive_sample_id(r1) == expected_sample


def test_derive_sample_id_returns_none_when_no_pattern_matches() -> None:
    assert derive_sample_id("A_R3.fastq.gz") is None


def test_derive_sample_id_allows_empty_prefix() -> None:
    # Ported from the bash script: `_L001_R1_001.fastq.gz` yields an empty
    # sample ID. (Pathological but should be observable to callers, not
    # silently wedged.)
    assert derive_sample_id("_L001_R1_001.fastq.gz") == ""


# ---------- discover() -----------------------------------------------------

def _make_fastq(folder: Path, name: str) -> Path:
    """Create an empty file with the given name; return its Path."""
    p = folder / name
    p.write_text("")
    return p


def test_discover_pairs_a_canonical_r1_with_its_r2(tmp_path: Path) -> None:
    r1 = _make_fastq(tmp_path, "Foo_L001_R1_001.fastq.gz")
    r2 = _make_fastq(tmp_path, "Foo_L001_R2_001.fastq.gz")

    result = discover([tmp_path])

    assert result == [Sample(sample_id="Foo", r1=r1, r2=r2)]


def test_discover_routes_single_end_when_no_pattern_matches(
    tmp_path: Path,
) -> None:
    f = _make_fastq(tmp_path, "loose_sample.fastq.gz")

    result = discover([tmp_path])

    assert result == [Sample(sample_id="loose_sample", r1=f, r2=None)]


def test_discover_routes_orphan_r1_as_single_end(tmp_path: Path) -> None:
    # An R1 whose R2 partner is missing is treated as single-end ([S21]),
    # per decision-4a. Sample ID still comes from extension-stripping
    # ([S21], not [S12]) — the pair convention does not apply when only
    # one mate is present.
    f = _make_fastq(tmp_path, "Foo_L001_R1_001.fastq.gz")

    result = discover([tmp_path])

    assert result == [Sample(sample_id="Foo_L001_R1_001", r1=f, r2=None)]


def test_discover_does_not_emit_r2_as_a_separate_sample(
    tmp_path: Path,
) -> None:
    # The R2 file must be claimed by its R1 partner, not surfaced as an
    # additional single-end sample.
    _make_fastq(tmp_path, "Foo_L001_R1_001.fastq.gz")
    _make_fastq(tmp_path, "Foo_L001_R2_001.fastq.gz")

    result = discover([tmp_path])

    assert len(result) == 1
    assert result[0].sample_id == "Foo"


def test_discover_mixes_paired_and_unpaired(tmp_path: Path) -> None:
    r1 = _make_fastq(tmp_path, "Paired_L001_R1_001.fastq.gz")
    r2 = _make_fastq(tmp_path, "Paired_L001_R2_001.fastq.gz")
    se = _make_fastq(tmp_path, "loose.fastq.gz")

    result = {s.sample_id: s for s in discover([tmp_path])}

    assert result == {
        "Paired": Sample(sample_id="Paired", r1=r1, r2=r2),
        "loose":  Sample(sample_id="loose",  r1=se, r2=None),
    }


def test_discover_walks_every_folder_in_the_list(tmp_path: Path) -> None:
    # [S10] — a list of folders is walked, every fastq file is collected.
    folder_a = tmp_path / "a"
    folder_b = tmp_path / "b"
    folder_a.mkdir()
    folder_b.mkdir()
    _make_fastq(folder_a, "From_A_L001_R1_001.fastq.gz")
    _make_fastq(folder_a, "From_A_L001_R2_001.fastq.gz")
    _make_fastq(folder_b, "from_b.fastq.gz")

    result = {s.sample_id for s in discover([folder_a, folder_b])}

    assert result == {"From_A", "from_b"}


def test_discover_accepts_every_extension(tmp_path: Path) -> None:
    extensions = (
        ".fastq", ".fq",
        ".fastq.gz", ".fq.gz",
        ".fastq.bz2", ".fq.bz2",
    )
    for ext in extensions:
        f = _make_fastq(tmp_path, f"sample{ext}")
        assert any(
            s.r1 == f and s.sample_id == "sample"
            for s in discover([tmp_path])
        ), f"missed extension {ext}"
        f.unlink()


def test_discover_ignores_non_fastq_files(tmp_path: Path) -> None:
    _make_fastq(tmp_path, "real.fastq.gz")
    (tmp_path / "notes.txt").write_text("ignore me")
    (tmp_path / "metadata.tsv").write_text("ignore me too")

    result = discover([tmp_path])

    assert [s.r1.name for s in result] == ["real.fastq.gz"]


# ---------- user --fastq_pattern override ----------------------------------

def test_discover_uses_user_pattern_to_pair_non_canonical_names(
    tmp_path: Path,
) -> None:
    # File names like `myrun_demoX_1.fq.gz` / `_2.fq.gz` do match
    # canonical row 7 already — pick a layout the canonical table
    # genuinely doesn't recognise to prove the override path.
    r1 = _make_fastq(tmp_path, "study42-mateA.fastq.gz")
    r2 = _make_fastq(tmp_path, "study42-mateB.fastq.gz")

    # Without an override, both files would look single-end:
    assert {s.sample_id for s in discover([tmp_path])} == {
        "study42-mateA", "study42-mateB",
    }

    # With the override, they pair up. The {1,2} token marks the
    # discriminator; the prefix `*` becomes the sample-ID capture.
    result = discover([tmp_path], extra_pattern="*-mate{A,B}.fastq.gz")

    assert result == [Sample(sample_id="study42", r1=r1, r2=r2)]


def test_user_pattern_without_star_uses_literal_prefix_as_sample_id(
    tmp_path: Path,
) -> None:
    # When the user pattern has no `*` in the prefix (e.g. the legacy
    # `params.fastq_pattern` of an end-to-end test that pins discovery
    # to a single named pair), the literal prefix minus trailing
    # separators is the sample ID — mirroring the historic
    # channel.fromFilePairs key-extraction behaviour.
    r1 = _make_fastq(tmp_path, "paired_merge_ok_1.fastq.gz")
    r2 = _make_fastq(tmp_path, "paired_merge_ok_2.fastq.gz")

    result = discover(
        [tmp_path], extra_pattern="paired_merge_ok_{1,2}.fastq.gz"
    )

    assert result == [Sample(sample_id="paired_merge_ok", r1=r1, r2=r2)]


def test_user_pattern_strips_legacy_leading_slash(tmp_path: Path) -> None:
    # `params.fastq_pattern` used to be concatenated onto fastq_folder
    # for fromFilePairs, so values commonly start with `/`. Treat the
    # input as a basename glob regardless.
    r1 = _make_fastq(tmp_path, "paired_merge_ok_1.fastq.gz")
    r2 = _make_fastq(tmp_path, "paired_merge_ok_2.fastq.gz")

    result = discover(
        [tmp_path], extra_pattern="/paired_merge_ok_{1,2}.fastq.gz"
    )

    assert result == [Sample(sample_id="paired_merge_ok", r1=r1, r2=r2)]


def test_user_pattern_takes_precedence_over_canonical_table(
    tmp_path: Path,
) -> None:
    # `A_R1.fastq.gz` / `A_R2.fastq.gz` would normally pair via canonical
    # row 6 with sample ID "A". A user pattern matching the same file
    # must win — this lets users override sample-ID derivation when the
    # canonical guess is wrong for their layout.
    r1 = _make_fastq(tmp_path, "A_R1.fastq.gz")
    r2 = _make_fastq(tmp_path, "A_R2.fastq.gz")

    result = discover([tmp_path], extra_pattern="*{R1,R2}.fastq.gz")

    # Override strips at the literal `R1`, so sample = "A_".
    assert result == [Sample(sample_id="A_", r1=r1, r2=r2)]


# ---------- sample-ID uniqueness ([S13]/[S14], D03) ------------------------
# COVERAGE: [S13], [S14]


def test_check_unique_sample_ids_passes_on_distinct_ids() -> None:
    samples = [
        Sample(sample_id="A", r1=Path("A_1.fastq.gz"), r2=None),
        Sample(sample_id="B", r1=Path("B_1.fastq.gz"), r2=None),
    ]
    # no exception
    check_unique_sample_ids(samples)


def test_check_unique_sample_ids_raises_when_two_samples_collide() -> None:
    samples = [
        Sample(
            sample_id="A",
            r1=Path("/run1/A_1.fastq.gz"),
            r2=Path("/run1/A_2.fastq.gz"),
        ),
        Sample(
            sample_id="A",
            r1=Path("/run2/A_R1.fastq.gz"),
            r2=Path("/run2/A_R2.fastq.gz"),
        ),
    ]
    with pytest.raises(DuplicateSampleIDError) as excinfo:
        check_unique_sample_ids(samples)
    msg = str(excinfo.value)
    # the colliding sample ID is named
    assert "A" in msg
    # every offending fastq path is listed (R1 paths at minimum)
    assert "/run1/A_1.fastq.gz" in msg
    assert "/run2/A_R1.fastq.gz" in msg


def test_check_unique_sample_ids_reports_every_collision_group() -> None:
    # two independent collisions must all surface in one error so the
    # user can resolve them in a single pass.
    samples = [
        Sample(sample_id="A", r1=Path("/run1/A_R1.fastq.gz"), r2=None),
        Sample(sample_id="A", r1=Path("/run2/A_R1.fastq.gz"), r2=None),
        Sample(sample_id="B", r1=Path("/run1/B.fastq.gz"), r2=None),
        Sample(sample_id="C", r1=Path("/run1/C_R1.fastq.gz"), r2=None),
        Sample(sample_id="C", r1=Path("/run2/C_R1.fastq.gz"), r2=None),
    ]
    with pytest.raises(DuplicateSampleIDError) as excinfo:
        check_unique_sample_ids(samples)
    msg = str(excinfo.value)
    assert "A" in msg and "C" in msg
    # the non-colliding sample (B) must not be listed
    assert "B.fastq.gz" not in msg


def test_check_unique_sample_ids_lists_both_r1_and_r2(tmp_path: Path) -> None:
    # When the colliding samples are paired-end, the error must surface
    # both R1 *and* R2 so the user knows which file pair to delete.
    samples = [
        Sample(
            sample_id="A",
            r1=Path("/run1/A_1.fastq.gz"),
            r2=Path("/run1/A_2.fastq.gz"),
        ),
        Sample(
            sample_id="A",
            r1=Path("/run2/A_1.fastq.gz"),
            r2=Path("/run2/A_2.fastq.gz"),
        ),
    ]
    with pytest.raises(DuplicateSampleIDError) as excinfo:
        check_unique_sample_ids(samples)
    msg = str(excinfo.value)
    for required in (
        "/run1/A_1.fastq.gz",
        "/run1/A_2.fastq.gz",
        "/run2/A_1.fastq.gz",
        "/run2/A_2.fastq.gz",
    ):
        assert required in msg, f"missing {required} in: {msg}"


def test_cli_exits_non_zero_on_duplicate_sample_ids(tmp_path: Path) -> None:
    # Two folders each providing a pair that resolves to sample ID "A".
    a = tmp_path / "run1"
    b = tmp_path / "run2"
    a.mkdir()
    b.mkdir()
    (a / "A_1.fastq.gz").write_bytes(b"")
    (a / "A_2.fastq.gz").write_bytes(b"")
    (b / "A_R1.fastq.gz").write_bytes(b"")
    (b / "A_R2.fastq.gz").write_bytes(b"")

    script = (
        Path(__file__).resolve().parents[2] / "bin" / "discover_fastq.py"
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
    # Both R1 paths must be listed
    assert str(a / "A_1.fastq.gz") in proc.stderr
    assert str(b / "A_R1.fastq.gz") in proc.stderr
