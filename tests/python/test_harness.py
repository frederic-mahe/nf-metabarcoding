"""Smoke tests for the pytest harness itself.

These tests have no `[Sxx]` coverage — they exist only to prove the
test harness is wired correctly so CI fails loudly if pytest cannot
discover or import anything. Delete this file once `bin/*.py` lands
and real tests replace it.

When real scripts arrive, the recommended patterns are:

1. **Importable module** (preferred for new code)::

       # bin/build_occurrence_table.py
       def collapse_clusters(samples): ...

       # tests/python/test_build_occurrence_table.py
       from build_occurrence_table import collapse_clusters

       def test_empty_sample_gets_zero_row():
           assert collapse_clusters({"sampleA": []}) == [
               {"sample": "sampleA", "abundance": 0},
           ]

2. **CLI subprocess** (fallback for legacy scripts)::

       import subprocess

       def test_cli_smoke(bin_dir):
           result = subprocess.run(
               [str(bin_dir / "old_script.py"), "--help"],
               capture_output=True, text=True, check=True,
           )
           assert "usage" in result.stdout.lower()
"""

from __future__ import annotations

import sys
from pathlib import Path


def test_repo_root_fixture_resolves(repo_root: Path) -> None:
    assert repo_root.is_dir()
    assert (repo_root / "main.nf").is_file()


def test_bin_dir_is_on_sys_path(bin_dir: Path) -> None:
    assert str(bin_dir) in sys.path
    assert bin_dir.is_dir()
