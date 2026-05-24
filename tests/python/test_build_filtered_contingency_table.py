"""Characterization tests for ``bin/build_filtered_contingency_table.py``.

This is the byte-exact golden-file gate that pins the legacy
script's current behaviour before any refactor. The fixture is
engineered to exercise every branch in a single run:

* cluster **aa01** — passes every filter; exercises:
  - ``;size=`` strip in ``stampa_parse`` taxonomy column (``Phy#lum1``)
  - duplicate sample rows in ``.distr`` (S1 appears twice for aa01 →
    summed)
  - a *non-seed* amplicon (aa02) whose distribution rows aggregate
    onto the seed
* cluster **ee01** — passes by ``abundance >= 3`` alone (``spread=1``)
* cluster **ff01** — passes by ``spread >= 2`` alone (``abundance < 3``)
* cluster **bb01** — filtered by ``chimera_status == "Y"``
* cluster **cc01** — filtered by ``ee/length > 0.0002``
* cluster **dd01** — filtered by ``abundance < 3 and spread < 2``
* cluster **1101** — filtered by ``chimera_status == "NA"`` because
  its uchime row has fewer than 18 columns (the script's second
  ``IndexError`` branch)
* an orphan **partial uchime line** at the top of ``.uchime`` to
  exercise the first ``IndexError`` branch (no ``OTU[1]``)
* an orphan **distribution row** (``9999``) for a sample ``s5`` that
  no surviving cluster references — ``s5`` must still appear as a
  zero column in the output
* sample **s4** lives only inside the filtered-out ``bb01`` row and
  must also appear as a zero column

If a future refactor changes the byte-exact output, this test
breaks; if the change is intentional, regenerate ``expected_table.tsv``
in the same commit and document the reason.
"""

# COVERAGE: [S35]

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


SCRIPT = (
    Path(__file__).resolve().parents[2]
    / "bin"
    / "build_filtered_contingency_table.py"
)
FIXTURE = (
    Path(__file__).resolve().parents[1]
    / "data"
    / "build_filtered_contingency_table"
)


def _run(extra: list[str] | None = None) -> subprocess.CompletedProcess[str]:
    cmd = [
        sys.executable,
        str(SCRIPT),
        "--representatives", str(FIXTURE / "representatives"),
        "--stats",           str(FIXTURE / "stats"),
        "--swarms",          str(FIXTURE / "swarms"),
        "--chimera",         str(FIXTURE / "uchime"),
        "--quality",         str(FIXTURE / "quality"),
        "--assignments",     str(FIXTURE / "assignments"),
        "--distribution",    str(FIXTURE / "distribution"),
    ]
    if extra:
        cmd += extra
    return subprocess.run(cmd, capture_output=True, text=True, check=True)


def test_golden_table() -> None:
    expected = (FIXTURE / "expected_table.tsv").read_text()
    assert _run().stdout == expected, "byte-exact mismatch"


def test_samples_flag_adds_zero_columns_for_empty_samples() -> None:
    """[S09]: ``--samples`` injects sample columns for inputs that
    contribute zero rows to the distribution file (entirely-empty
    samples). The fixture's distribution lists s1..s5; adding sNew
    must surface sNew as a zero column on every output row.

    Back-compat: samples already in the distribution stay unaffected
    (union with the --samples set, no duplicates).
    """
    result = _run(["--samples", "sNew,s2,sExtra"])
    lines = result.stdout.splitlines()
    header = lines[0].split("\t")

    # Sample columns occupy cols 13 onward (after the 13 metadata cols).
    sample_cols = header[13:]
    assert sample_cols == ["s1", "s2", "s3", "s4", "s5", "sExtra", "sNew"], (
        f"expected union of distribution samples + --samples, sorted; "
        f"got: {sample_cols}"
    )

    # Every data row carries zero for the newly-injected samples.
    new_indices = {name: 13 + sample_cols.index(name)
                   for name in ("sNew", "sExtra")}
    for row in lines[1:]:
        cols = row.split("\t")
        for name, idx in new_indices.items():
            assert cols[idx] == "0", (
                f"injected sample {name} must be 0 on every row, got "
                f"'{cols[idx]}' in row: {row}"
            )
