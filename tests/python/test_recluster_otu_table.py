"""Characterization tests for ``bin/recluster_otu_table.py``.

Byte-exact golden-file gate for the optional post-mumu re-clustering
merge ([S104]/[S105], D20). ``recluster_otu_table.py`` is a cousin of
``merge_substring_otus.py`` kept separate so [S39]'s golden output
stays byte-exact; the differences it pins here are:

* it keys on the ``amplicon`` column and strips the ``;size=N;``
  annotation vsearch ``--cluster_size --sizeout`` writes onto the
  ``.uc`` labels;
* ``cloud`` is left at the post-mumu ``"NA"`` string (no ``+1`` per
  merged pupil — that is the substring-merge quirk, not this one);
* surviving centroid rows are renumbered ``1..N`` in input order
  (D-d), so no downstream bash sort is needed.

The main fixture exercises: a pass-through singleton (ampF), a centroid
with one member (ampB←ampD), a centroid with two members
(ampA←ampC,ampE), ``cloud`` staying ``NA``, ``spread`` recomputation,
contiguous renumbering (6 OTUs → 3), and read-count conservation.

A separate test pins the overlap abort (an OTU that is both centroid
and member), and an inline test pins that a member row appearing
*before* its centroid in the table still folds correctly.

The match file format mirrors vsearch's ``--uc`` ``^H`` lines: the
last two tab-separated fields are ``<member>\\t<centroid>``, each a
``<amplicon>;size=N;`` label.
"""

# COVERAGE: [S104]

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


SCRIPT = (
    Path(__file__).resolve().parents[2]
    / "bin"
    / "recluster_otu_table.py"
)
FIXTURE = (
    Path(__file__).resolve().parents[1]
    / "data"
    / "recluster_otu_table"
)


def test_golden_output(tmp_path: Path) -> None:
    out = tmp_path / "reclustered_table"
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


def test_read_count_conserved(tmp_path: Path) -> None:
    out = tmp_path / "reclustered_table"
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

    def total_sum(path: Path) -> int:
        lines = path.read_text().splitlines()[1:]
        return sum(int(line.split("\t")[1]) for line in lines)

    assert total_sum(out) == total_sum(FIXTURE / "otu_table")


def test_centroid_member_overlap_aborts(tmp_path: Path) -> None:
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


def test_member_before_centroid_folds(tmp_path: Path) -> None:
    """A member row that precedes its centroid still folds in.

    vsearch orders clusters by centroid abundance, but the post-mumu
    table is ordered by OTU number ([S44]), so a low-abundance member
    can sit above its centroid. The merge must not depend on row order.
    """
    header = (
        "OTU\ttotal\tcloud\tamplicon\tlength\tabundance\tchimera\t"
        "spread\tquality\tsequence\tidentity\ttaxonomy\treferences\ts1\n"
    )
    # OTU 1 (ampLo) is a member of OTU 2 (ampHi); the member appears first.
    table = tmp_path / "table"
    table.write_text(
        header
        + "1\t3\tNA\tampLo\t50\t3\tN\t1\t0.001\tAC\t0.0\tNA\tNA\t3\n"
        + "2\t9\tNA\tampHi\t50\t9\tN\t1\t0.001\tGT\t0.0\tNA\tNA\t9\n"
    )
    matches = tmp_path / "matches"
    matches.write_text(
        "H\t*\t*\t*\t*\t*\t*\t*\tampLo;size=3;\tampHi;size=9;\n"
    )
    out = tmp_path / "out"
    subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "-t", str(table),
            "-m", str(matches),
            "-o", str(out),
        ],
        check=True,
    )
    rows = out.read_text().splitlines()
    assert len(rows) == 2, f"expected header + 1 centroid, got: {rows}"
    fields = rows[1].split("\t")
    assert fields[0] == "1", "surviving centroid should be renumbered to 1"
    assert fields[3] == "ampHi", "centroid metadata should win"
    assert fields[1] == "12", "total should be 9 + 3"
    assert fields[-1] == "12", "sample column should be 9 + 3"
