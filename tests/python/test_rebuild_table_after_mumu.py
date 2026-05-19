"""Characterization tests for ``bin/rebuild_table_after_mumu.py``.

Byte-exact golden-file gate. The fixture exercises:

* the header passthrough (line 0 of the old table → stdout verbatim)
* per-amplicon metadata join (``length``/``abundance``/``quality``/
  ``sequence``/``identity``/``taxonomy``/``references`` from the old
  table) keyed on the amplicon ID at column 3
* per-row recomputation: new ``total`` = sum of mumu sample columns,
  new ``spread`` = count of non-zero sample columns
* the ``cloud`` and ``chimera`` placeholders (``"NA"`` and ``"N"``
  constants — mumu drops the cloud info)
* a zero-abundance row surviving the join (``dd04``: total=0,
  spread=0); the size=0 → 1 hotfix is applied downstream by the
  bash awk wrapper, not by this script.

OTU numbering restarts at 1 and increments with each emitted mumu
row, regardless of the input OTU numbers (which would have come
from the pre-mumu table).
"""

# COVERAGE: [S44]

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


SCRIPT = (
    Path(__file__).resolve().parents[2]
    / "bin"
    / "rebuild_table_after_mumu.py"
)
FIXTURE = (
    Path(__file__).resolve().parents[1]
    / "data"
    / "rebuild_table_after_mumu"
)


def test_golden_output() -> None:
    cmd = [
        sys.executable,
        str(SCRIPT),
        "--mumu_table", str(FIXTURE / "mumu_table"),
        "--old_table",  str(FIXTURE / "old_table"),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    expected = (FIXTURE / "expected_output").read_text()
    assert result.stdout == expected, (
        "byte-exact mismatch.\n"
        f"Got:\n{result.stdout!r}\n"
        f"Expected:\n{expected!r}"
    )
