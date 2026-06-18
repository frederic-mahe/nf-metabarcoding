#!/usr/bin/env python3
"""Collect external-tool version strings into a YAML mapping ([S68]).

Reproducibility: every nf-metabarcoding run records the versions of the
external tools it used (``vsearch``, ``swarm``, ``cutadapt``, ``mumu``)
plus the Python interpreter into
``<results_folder>/pipeline_info/software_versions.yml``.

Each input line is ``<name>\\t<raw version output>`` — the raw text a
tool prints for ``--version``. The version token is extracted with a
single permissive regex (the first ``MAJOR.MINOR[.PATCH...]`` run), so
the wildly different shapes the tools print (``Swarm 3.1.6``,
``vsearch v2.31.0_linux_x86_64, ...``, a bare ``5.2``) all reduce to a
clean version. A tool that prints nothing recognisable (e.g. missing
from PATH, ``command not found``) is recorded as ``n/a`` rather than
dropped, so a broken environment is visible in the report instead of
silently absent.

Reads the lines from stdin (or from files passed as arguments) and
writes the sorted YAML mapping to stdout.
"""

from __future__ import annotations

import re
import sys
from typing import Iterable

_VERSION_RE = re.compile(r"\d+\.\d+(?:\.\d+)*")


def extract_version(raw: str) -> str:
    """Return the first ``MAJOR.MINOR[.PATCH...]`` token, or ``n/a``."""
    match = _VERSION_RE.search(raw)
    return match.group(0) if match else "n/a"


def collect(lines: Iterable[str]) -> dict[str, str]:
    """Map each ``name<TAB>raw`` line to ``name -> extracted version``.

    Blank lines and lines with no name are skipped; later duplicates
    of the same name win.
    """
    versions: dict[str, str] = {}
    for line in lines:
        if not line.strip():
            continue
        name, _, raw = line.rstrip("\n").partition("\t")
        name = name.strip()
        if not name:
            continue
        versions[name] = extract_version(raw)
    return versions


def to_yaml(versions: dict[str, str]) -> str:
    """Render *versions* as a deterministic, sorted YAML mapping."""
    return "".join(
        f"{name}: {versions[name]}\n" for name in sorted(versions)
    )


def main(argv: list[str] | None = None) -> int:
    args = sys.argv[1:] if argv is None else argv
    if args:
        lines: list[str] = []
        for path in args:
            with open(path, encoding="utf-8") as handle:
                lines.extend(handle.readlines())
    else:
        lines = sys.stdin.readlines()
    sys.stdout.write(to_yaml(collect(lines)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
