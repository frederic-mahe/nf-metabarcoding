"""Every ``bin/*.py`` helper fails early on a too-old Python ([S109]).

Without an engine profile, tasks run the host's ``python3`` (Python
3.6 on the Genotoul compute nodes). The helpers must then stop with a
readable message rather than a ``SyntaxError`` deep in the file:

* no ``from __future__ import annotations`` (a compile-time error on
  Python < 3.7) and a source that parses under the 3.6 grammar, so the
  version check is always reached;
* with ``sys.version_info`` patched to 3.6, running the script exits 1
  and names the script and the required version on stderr.
"""

# COVERAGE: [S109]

from __future__ import annotations

import ast
import subprocess
import sys
from pathlib import Path

import pytest

BIN_DIR = Path(__file__).resolve().parents[2] / "bin"
SCRIPTS = sorted(BIN_DIR.glob("*.py"))

# Run the script as __main__ under a faked Python 3.6.
FAKE_OLD_PYTHON = (
    "import runpy, sys\n"
    "sys.version_info = (3, 6, 8, 'final', 0)\n"
    "sys.version = '3.6.8 (default)'\n"
    "runpy.run_path(sys.argv[1], run_name='__main__')\n"
)


def test_scripts_found() -> None:
    assert SCRIPTS


@pytest.mark.parametrize("script", SCRIPTS, ids=lambda p: p.name)
def test_no_future_annotations(script: Path) -> None:
    assert "from __future__ import annotations" not in script.read_text()


@pytest.mark.parametrize("script", SCRIPTS, ids=lambda p: p.name)
def test_parses_under_python_3_6_grammar(script: Path) -> None:
    ast.parse(script.read_text(), feature_version=(3, 6))


@pytest.mark.parametrize("script", SCRIPTS, ids=lambda p: p.name)
def test_exits_early_on_old_python(script: Path) -> None:
    result = subprocess.run(
        [sys.executable, "-c", FAKE_OLD_PYTHON, str(script)],
        capture_output=True,
        text=True,
        cwd=BIN_DIR,
    )
    assert result.returncode == 1
    assert "Python >= 3.10 required" in result.stderr
    assert "found 3.6.8" in result.stderr
    assert script.name in result.stderr
    assert "Traceback" not in result.stderr
