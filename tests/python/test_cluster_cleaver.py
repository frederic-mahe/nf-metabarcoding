"""Characterization tests for ``bin/cluster_cleaver.py``.

Golden-file tests pinned to the byte-exact output of the legacy
``tmp/OTU_cleaver.py`` on a hand-crafted fixture. The fixture is
engineered to exercise the script's branches in a single run:

* cluster A — splits along a sub-seed (basic ``[S22]`` path)
* cluster D — splits with a *sort tie* in the root sub-cluster
  (``dd01`` and ``dd04`` both at abundance 20; tie broken by name)
* cluster B — no sub-seed survives → cluster is NOT re-emitted
* cluster C — candidate ``cc02`` is below the per-sample threshold
  (1 out of 21 samples; threshold = 0.05 * 21 = 1.05) → cluster is
  NOT re-emitted
* cluster E — trivial split (sub-seed at the only struct line)
* cluster F — splits with a *seed change after sort*: the surviving
  sub-cluster's original key ``ff02`` (abundance 5) loses to ``ff03``
  (abundance 8) once amplicon abundances are filled in
* cluster Y — singleton (``cloud=1``); exercises the stats-parse
  filter that skips clusters too small to host a sub-seed
* cluster Z — trailing splitter present so cluster F's father-son
  walk completes before the ``number_of_seeds==0`` early break fires
* cluster W — final dummy whose only role is to trigger the cluster
  transition that appends Z and then fires the early-break exit

Two test cases cover both filename-driven output suffixes: the
``_1f`` (fastidious) path that matches ``swarm --fastidious`` runs,
and the non-fastidious path used when neither the swarms nor struct
file name contains ``_1f.``.
"""

# COVERAGE: [S22]

from __future__ import annotations

import subprocess
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[2] / "bin" / "cluster_cleaver.py"


# ---------- fixture inputs -------------------------------------------------

# Amplicon names use 2-hex-char prefixes (``aa``..``ff``, ``00``..``ff``)
# because the legacy script hashes each name into 256 buckets via
# ``int(name[0:2], 16)``.
INPUT_FASTA = """\
>aa01;size=30
AAAAAAAAAA
>ff01;size=30
TGAGAGAGTC
>dd01;size=20
GGGGGGGGGG
>dd04;size=20
GGGGGCGGGG
>bb01;size=15
CCCCCCCCCC
>cc01;size=12
TTTTTTTTTT
>0001;size=10
TAGCTAGCTA
>ee01;size=8
GAGAGAGAGA
>ff03;size=8
TGAGAGAGTT
>0a01;size=6
AGTCAGTCAG
>aa02;size=5
ATAAAAAAAA
>dd02;size=5
GTGGGGGGGG
>ff02;size=5
TTAGAGAGTC
>ee02;size=4
GTGAGAGAGA
>0a02;size=4
AGTCAGTCAA
>bb02;size=3
CTCCCCCCCC
>aa03;size=2
ATAAAAAAAT
>cc02;size=2
TATTTTTTTT
>1101;size=2
ACTGACTGAC
>1102;size=2
ACTGACTGAA
>cc03;size=1
TATTTTTTTA
"""

# One row per cluster (cluster id == row index, 1-based). Columns:
# cloud, mass, seed, seed_abundance, singletons. Row 7 (cluster Y)
# has cloud=1, exercising the stats-parse filter skip.
INPUT_STATS = """\
3\t37\taa01\t30\t0
3\t45\tdd01\t20\t0
2\t18\tbb01\t15\t0
3\t15\tcc01\t12\t1
2\t12\tee01\t8\t0
3\t43\tff01\t30\t0
1\t10\t0001\t10\t0
2\t4\t1101\t2\t0
2\t10\t0a01\t6\t0
"""

# Father-son tree. Columns: father, son, diffs, cluster_id, steps.
# Cluster 7 (singleton Y) has no struct lines. Each splitting cluster
# is shaped so its surviving sub-seed appears at the LAST struct line
# of the cluster, so the next cluster's transition appends it before
# the ``number_of_seeds==0`` early break fires.
INPUT_STRUCT = """\
aa01\taa03\t1\t1\t1
aa03\taa02\t1\t1\t2
dd01\tdd04\t1\t2\t1
dd01\tdd02\t1\t2\t2
bb01\tbb02\t1\t3\t1
cc01\tcc02\t1\t4\t1
cc02\tcc03\t1\t4\t2
ee01\tee02\t1\t5\t1
ff01\tff02\t1\t6\t1
ff02\tff03\t1\t6\t2
1101\t1102\t1\t8\t1
0a01\t0a02\t1\t9\t1
"""

# One line per cluster (cluster id == line index, 1-based). Format:
# ``<amplicon>;size=<n>`` space-separated.
INPUT_SWARMS = """\
aa01;size=30 aa02;size=5 aa03;size=2
dd01;size=20 dd04;size=20 dd02;size=5
bb01;size=15 bb02;size=3
cc01;size=12 cc02;size=2 cc03;size=1
ee01;size=8 ee02;size=4
ff01;size=30 ff02;size=5 ff03;size=8
0001;size=10
1101;size=2 1102;size=2
0a01;size=6 0a02;size=4
"""


# Per-sample cluster seeds. With 21 samples the threshold is
# 0.05 * 21 = 1.05, so a seed must appear in at least 2 samples to
# survive. ``aa02`` covers the "kept w/ count > 1" branch, ``cc02``
# covers the "dropped below threshold" branch, and the others
# (``dd02``, ``ee02``, ``ff02``, ``1102``) appear in every sample
# so that cluster F's walk can still find a member after its own
# sub-seed has been found.
def _build_per_sample_stats() -> str:
    sample_to_seeds: dict[str, list[str]] = {}
    for i in range(1, 6):
        sample_to_seeds[f"sample_{i:02d}"] = [
            "aa02", "dd02", "ee02", "ff02", "1102",
        ]
    sample_to_seeds["sample_06"] = ["cc02", "dd02", "ee02", "ff02", "1102"]
    for i in range(7, 22):
        sample_to_seeds[f"sample_{i:02d}"] = [
            "dd02", "ee02", "ff02", "1102",
        ]

    # cloud/mass/abundance values below are arbitrary — the legacy
    # script only reads columns 0 (sample) and 3 (seed).
    seed_to_stats: dict[str, tuple[str, str, str]] = {
        "aa02": ("2", "8", "5"),
        "cc02": ("2", "5", "2"),
        "dd02": ("3", "30", "20"),
        "ee02": ("2", "12", "4"),
        "ff02": ("3", "43", "5"),
        "1102": ("2", "4", "2"),
    }

    lines: list[str] = []
    for sample in sorted(sample_to_seeds):
        for seed in sample_to_seeds[sample]:
            cloud, mass, abundance = seed_to_stats[seed]
            lines.append(
                "\t".join((sample, cloud, mass, seed, abundance))
            )
    return "\n".join(lines) + "\n"


PER_SAMPLE_STATS = _build_per_sample_stats()


# ---------- expected outputs ----------------------------------------------

# Sub-cluster stats — sorted first by seed name (stable, asc), then
# by (total_abundance, num_uniques) desc. Cluster B and C are absent
# because no surviving sub-seed mapped into them. The ``2 32`` row
# is cluster A's root sub-cluster (aa01 + aa03 after the orphan-style
# father-son walk); the ``2 13 ff03`` row demonstrates the
# seed-change-after-sort path.
EXPECTED_STATS2 = """\
2\t40\tdd01\t20\t0\t0\t0
2\t32\taa01\t30\t0\t0\t0
1\t30\tff01\t30\t0\t0\t0
2\t13\tff03\t8\t0\t0\t0
1\t8\tee01\t8\t0\t0\t0
1\t5\taa02\t5\t0\t0\t0
1\t5\tdd02\t5\t0\t0\t0
1\t4\tee02\t4\t0\t0\t0
1\t2\t1101\t2\t0\t0\t0
1\t2\t1102\t2\t0\t0\t0
"""

# Sub-cluster swarms, emitted in ``new_clusters`` traversal order
# (cluster A, D, E, F, Z). Within each sub-cluster the amplicons are
# sorted by ``(-abundance, name)``: the ``dd01;size=20 dd04;size=20``
# line proves that the tie is broken by name; the
# ``ff03;size=8 ff02;size=5`` line proves the sub-cluster key swapped
# from ff02 to ff03 after the sort.
EXPECTED_SWARMS2 = """\
aa01;size=30 aa03;size=2
aa02;size=5
dd01;size=20 dd04;size=20
dd02;size=5
ee01;size=8
ee02;size=4
ff01;size=30
ff03;size=8 ff02;size=5
1101;size=2
1102;size=2
"""

# Sub-cluster representative FASTA — emitted in EXPECTED_STATS2 row
# order, with ``size=`` set to the sub-cluster's *total* abundance
# (not the seed's individual abundance).
EXPECTED_FAS2 = """\
>dd01;size=40
GGGGGGGGGG
>aa01;size=32
AAAAAAAAAA
>ff01;size=30
TGAGAGAGTC
>ff03;size=13
TGAGAGAGTT
>ee01;size=8
GAGAGAGAGA
>aa02;size=5
ATAAAAAAAA
>dd02;size=5
GTGGGGGGGG
>ee02;size=4
GTGAGAGAGA
>1101;size=2
ACTGACTGAC
>1102;size=2
ACTGACTGAA
"""


# ---------- helpers --------------------------------------------------------

def _write_inputs(
    folder: Path, *, swarm_parameters: str
) -> dict[str, Path]:
    """Materialize the embedded fixtures inside ``folder``.

    ``swarm_parameters`` is either ``"1f"`` (fastidious, mirrors the
    legacy pipeline) or ``"1"`` — it only shapes the input filenames
    so the script's ``"_1f." in swarms_file`` detection picks the
    right output FASTA suffix.
    """
    stem = f"cleaver_{swarm_parameters}" if swarm_parameters else "cleaver"
    paths = {
        "fasta": folder / "cleaver.fasta",
        "stats": folder / f"{stem}.stats",
        "struct": folder / f"{stem}.struct",
        "swarms": folder / f"{stem}.swarms",
        "per_sample": folder / "cleaver_per_sample.stats",
    }
    paths["fasta"].write_text(INPUT_FASTA)
    paths["stats"].write_text(INPUT_STATS)
    paths["struct"].write_text(INPUT_STRUCT)
    paths["swarms"].write_text(INPUT_SWARMS)
    paths["per_sample"].write_text(PER_SAMPLE_STATS)
    return paths


def _run_cleaver(
    inputs: dict[str, Path],
    *,
    cwd: Path,
    fastidious: bool = True,
    percentage: float | None = None,
) -> None:
    """Invoke ``bin/cluster_cleaver.py``.

    The script defaults to ``--fastidious`` (the legacy pipeline always
    runs swarm with ``--fastidious``). Pass ``fastidious=False`` to
    exercise the ``_1_`` filename branch. ``percentage`` overrides the
    sub-seed presence threshold (default 0.05; see [S22]).
    """
    cmd = [
        "python3",
        str(SCRIPT),
        "--global_stats", str(inputs["stats"]),
        "--per_sample_stats", str(inputs["per_sample"]),
        "--fasta", str(inputs["fasta"]),
        "--struct", str(inputs["struct"]),
        "--swarms", str(inputs["swarms"]),
    ]
    if not fastidious:
        cmd.append("--no-fastidious")
    if percentage is not None:
        cmd += ["--percentage", str(percentage)]
    subprocess.run(cmd, cwd=cwd, check=True)


# ---------- tests ----------------------------------------------------------

def test_cleaver_fastidious_default(tmp_path: Path) -> None:
    """End-to-end golden test: default --fastidious yields '_1f_' suffix."""
    inputs = _write_inputs(tmp_path, swarm_parameters="1f")
    _run_cleaver(inputs, cwd=tmp_path)

    assert (tmp_path / "cleaver_1f.stats2").read_text() == EXPECTED_STATS2
    assert (tmp_path / "cleaver_1f.swarms2").read_text() == EXPECTED_SWARMS2
    assert (
        tmp_path / "cleaver_1f_representatives.fas2"
    ).read_text() == EXPECTED_FAS2


def test_cleaver_no_fastidious_flag(tmp_path: Path) -> None:
    """Passing --no-fastidious selects the '_1_' representative suffix.

    The choice is now driven by the explicit flag (default True), not
    by ``_1f.`` appearing in the input filenames — so the stats and
    swarms files retain their original names even when the flag flips.
    """
    inputs = _write_inputs(tmp_path, swarm_parameters="1f")
    _run_cleaver(inputs, cwd=tmp_path, fastidious=False)

    assert (tmp_path / "cleaver_1f.stats2").read_text() == EXPECTED_STATS2
    assert (tmp_path / "cleaver_1f.swarms2").read_text() == EXPECTED_SWARMS2
    assert (
        tmp_path / "cleaver_1_representatives.fas2"
    ).read_text() == EXPECTED_FAS2
    # The fastidious-suffix variant must NOT exist when --no-fastidious
    # is passed — guards the flag-driven branch in ``main()``.
    assert not (tmp_path / "cleaver_1f_representatives.fas2").exists()


def test_cleaver_no_subseeds_emits_empty_outputs(tmp_path: Path) -> None:
    """Fixtures where no cluster has a sub-seed must still produce
    output files (empty), not crash.

    Mirrors what happens when the Part B pipeline runs on small data:
    if the per-sample stats never surface a sub-seed above the
    threshold, ``new_stats`` ends up empty and the cleaver historically
    crashed inside ``fasta_parse`` with ``ValueError: min() iterable
    argument is empty``. Empty outputs are a legitimate outcome
    (downstream concatenation via ``cat <out>{,2}`` tolerates empty
    files), so the script must produce them rather than abort.
    """
    fasta = tmp_path / "trivial.fasta"
    stats = tmp_path / "trivial_1f.stats"
    struct = tmp_path / "trivial_1f.struct"
    swarms = tmp_path / "trivial_1f.swarms"
    per_sample = tmp_path / "trivial_per_sample.stats"

    # Single global cluster, single per-sample seed, but the sample's
    # seed is the same amplicon as the global seed → no sub-seed
    # candidate, no cleaving.
    fasta.write_text(">aa01;size=10\nACGTACGT\n")
    stats.write_text("1\t10\taa01\t10\t0\t0\t0\n")
    struct.write_text("")
    swarms.write_text("aa01;size=10\n")
    # one sample row, with the global seed as its own per-sample seed
    per_sample.write_text("sample1\t1\t10\taa01\t10\t0\t0\t0\n")

    subprocess.run(
        [
            "python3",
            str(SCRIPT),
            "--global_stats", str(stats),
            "--per_sample_stats", str(per_sample),
            "--fasta", str(fasta),
            "--struct", str(struct),
            "--swarms", str(swarms),
        ],
        cwd=tmp_path,
        check=True,
    )

    # All three legacy `<input>2` files exist and are empty.
    out_stats = tmp_path / "trivial_1f.stats2"
    out_swarms = tmp_path / "trivial_1f.swarms2"
    out_fasta = tmp_path / "trivial_1f_representatives.fas2"
    assert out_stats.exists()
    assert out_swarms.exists()
    assert out_fasta.exists()
    assert out_stats.read_text() == ""
    assert out_swarms.read_text() == ""
    assert out_fasta.read_text() == ""


def test_cleaver_explicit_out_paths(tmp_path: Path) -> None:
    """``--out-stats / --out-swarms / --out-fasta`` redirect outputs."""
    inputs = _write_inputs(tmp_path, swarm_parameters="1f")
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    out_stats = out_dir / "custom.stats"
    out_swarms = out_dir / "custom.swarms"
    out_fasta = out_dir / "custom.fasta"

    subprocess.run(
        [
            "python3",
            str(SCRIPT),
            "--global_stats", str(inputs["stats"]),
            "--per_sample_stats", str(inputs["per_sample"]),
            "--fasta", str(inputs["fasta"]),
            "--struct", str(inputs["struct"]),
            "--swarms", str(inputs["swarms"]),
            "--out-stats", str(out_stats),
            "--out-swarms", str(out_swarms),
            "--out-fasta", str(out_fasta),
        ],
        cwd=tmp_path,
        check=True,
    )

    assert out_stats.read_text() == EXPECTED_STATS2
    assert out_swarms.read_text() == EXPECTED_SWARMS2
    assert out_fasta.read_text() == EXPECTED_FAS2
    # Legacy '<input>2'-suffixed files must NOT be created when the
    # caller passed explicit output paths.
    assert not (tmp_path / "cleaver_1f.stats2").exists()
    assert not (tmp_path / "cleaver_1f.swarms2").exists()
    assert not (tmp_path / "cleaver_1f_representatives.fas2").exists()


def test_cleaver_percentage_raises_threshold(tmp_path: Path) -> None:
    """[S22]: ``--percentage`` controls the sub-seed presence threshold.

    The PER_SAMPLE_STATS fixture spans 21 samples. At the default
    ``--percentage 0.05`` the threshold is 1.05, so any sub-seed
    appearing in ≥2 samples survives. Five clusters end up cleaved
    (A, D, E, F, 11 — see ``EXPECTED_SWARMS2``).

    Sub-seed presence counts in the fixture:
        aa02:  5 samples   (samples 01–05)
        cc02:  1 sample    (sample 06; dropped at every threshold)
        dd02:  21 samples  (all samples)
        ee02:  21 samples  (all samples)
        ff02:  21 samples  (all samples)
        1102:  21 samples  (all samples)

    Raising ``--percentage`` to 0.30 lifts the threshold to 6.3 — only
    ``aa02`` falls below it. Cluster A's only sub-seed is dropped, so
    cluster A is NOT re-emitted; clusters D / E / F / 11 still split.
    Four cleaved clusters × 2 sub-clusters each = 8 swarms2 lines
    (down from the default's 10).

    A regression that hard-coded the threshold (or ignored the flag)
    would still emit cluster A's two sub-cluster lines and fail both
    the count assertion and the per-key check.
    """
    inputs = _write_inputs(tmp_path, swarm_parameters="1f")
    _run_cleaver(inputs, cwd=tmp_path, percentage=0.30)

    swarms2 = (tmp_path / "cleaver_1f.swarms2").read_text().splitlines()
    assert len(swarms2) == 8, (
        f"expected 8 sub-cluster lines under --percentage 0.30, got "
        f"{len(swarms2)}:\n" + "\n".join(swarms2)
    )

    # aa02 is the only sub-seed dropped at this threshold. It must not
    # lead a sub-cluster line (i.e. it never appears as the surviving
    # sub-cluster's key) and the parent aa01 cluster must therefore be
    # absent from the output entirely.
    leading_keys = {line.split(";", 1)[0] for line in swarms2}
    assert "aa02" not in leading_keys, (
        f"aa02 should be filtered out at --percentage 0.30 "
        f"(5/21 samples = 23.8%), got sub-cluster keys: {leading_keys}"
    )
    assert "aa01" not in leading_keys, (
        f"cluster A has no surviving sub-seed at --percentage 0.30; its "
        f"master sub-cluster (keyed on aa01) must not be emitted, got "
        f"sub-cluster keys: {leading_keys}"
    )

    # The 21/21 sub-seeds still drive their splits.
    for kept_seed in ("dd02", "ee02", "ff03", "1102"):
        assert kept_seed in leading_keys, (
            f"{kept_seed} should still drive a split at --percentage "
            f"0.30; got sub-cluster keys: {leading_keys}"
        )
