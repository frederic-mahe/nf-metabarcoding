"""Unit tests for ``bin/collect_versions.py`` ([S68]).

The helper turns each external tool's raw ``--version`` output into a
clean, sorted YAML mapping for
``pipeline_info/software_versions.yml``. The extraction must cope with
the very different shapes the four external tools (and the Python
interpreter) print, and must record a tool that is missing from the
environment as ``n/a`` rather than dropping it — a silent omission
would hide a broken environment.
"""

# COVERAGE: [S68]

from __future__ import annotations

import subprocess
import sys

import pytest

from collect_versions import collect, extract_version, to_yaml


# Real first-line `--version` output captured from the pinned tools.
REAL_OUTPUTS = {
    "vsearch": "vsearch v2.31.0_linux_x86_64, 125.3GB RAM, 24 cores",
    "swarm": "Swarm 3.1.6",
    "cutadapt": "5.2",
    "mumu": "mumu 1.1.3",
    "python": "Python 3.12.3",
}
EXPECTED = {
    "vsearch": "2.31.0",
    "swarm": "3.1.6",
    "cutadapt": "5.2",
    "mumu": "1.1.3",
    "python": "3.12.3",
}


@pytest.mark.parametrize("tool", sorted(REAL_OUTPUTS))
def test_extract_version_per_tool(tool: str) -> None:
    assert extract_version(REAL_OUTPUTS[tool]) == EXPECTED[tool]


def test_extract_ignores_trailing_build_metadata() -> None:
    # vsearch appends `_linux_x86_64` and RAM/core counts; only the
    # MAJOR.MINOR.PATCH token is kept.
    raw = "vsearch v2.31.0_linux_x86_64, 125.3GB RAM"
    assert extract_version(raw) == "2.31.0"


def test_missing_tool_is_recorded_as_na() -> None:
    assert extract_version("") == "n/a"
    assert extract_version("vsearch: command not found") == "n/a"
    assert collect(["mumu\t"]) == {"mumu": "n/a"}


def test_collect_maps_name_to_version() -> None:
    lines = [f"{t}\t{REAL_OUTPUTS[t]}" for t in REAL_OUTPUTS]
    assert collect(lines) == EXPECTED


def test_collect_skips_blank_and_nameless_lines() -> None:
    assert collect(["", "   ", "\tno-name-here"]) == {}


def test_to_yaml_is_sorted_and_terminated() -> None:
    out = to_yaml({"swarm": "3.1.6", "cutadapt": "5.2"})
    assert out == "cutadapt: 5.2\nswarm: 3.1.6\n"


def test_cli_reads_stdin(bin_dir) -> None:
    script = bin_dir / "collect_versions.py"
    stdin = "".join(f"{t}\t{REAL_OUTPUTS[t]}\n" for t in REAL_OUTPUTS)
    proc = subprocess.run(
        [sys.executable, str(script)],
        input=stdin,
        capture_output=True,
        text=True,
        check=True,
    )
    assert proc.stdout == to_yaml(EXPECTED)
