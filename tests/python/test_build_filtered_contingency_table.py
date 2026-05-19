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


def test_golden_table(tmp_path: Path) -> None:
    cmd = [
        sys.executable,
        str(SCRIPT),
        "--representatives", str(FIXTURE / "representatives.fas"),
        "--stats",           str(FIXTURE / "stats"),
        "--swarms",          str(FIXTURE / "swarms"),
        "--chimera",         str(FIXTURE / "uchime"),
        "--quality",         str(FIXTURE / "quality"),
        "--assignments",     str(FIXTURE / "assignments"),
        "--distribution",    str(FIXTURE / "distribution"),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)

    expected = (FIXTURE / "expected_table.tsv").read_text()
    assert result.stdout == expected, (
        "byte-exact mismatch. Got:\n"
        f"{result.stdout!r}\n"
        f"Expected:\n{expected!r}"
    )
