"""Test configuration for pytest.

Adds the repository's ``bin/`` directory to ``sys.path`` so tests can
import helper modules (e.g. ``import build_occurrence_table``)
without packaging the workflow as a wheel. CLI-only scripts can still
be tested via :func:`subprocess.run` — see ``test_harness.py`` for
the recommended pattern.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest


REPO_ROOT = Path(
    subprocess.check_output(
        ["git", "rev-parse", "--show-toplevel"], text=True
    ).strip()
)
BIN_DIR = REPO_ROOT / "bin"

if str(BIN_DIR) not in sys.path:
    sys.path.insert(0, str(BIN_DIR))


@pytest.fixture(scope="session")
def repo_root() -> Path:
    """Absolute path to the repository root."""
    return REPO_ROOT


@pytest.fixture(scope="session")
def bin_dir() -> Path:
    """Absolute path to the ``bin/`` directory."""
    return BIN_DIR
