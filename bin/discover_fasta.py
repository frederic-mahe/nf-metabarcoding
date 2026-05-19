#!/usr/bin/env python3
"""Discover ``.fas`` files for Part B's fasta channel.

Used as:

* an importable module — call :func:`discover` or
  :func:`check_unique_sample_ids` directly;
* a CLI — ``discover_fasta.py FOLDER [FOLDER ...]`` emits TSV
  ``sample_id\\tfasta_path`` to stdout.

Rules (see ``[S27]`` in SPECIFICATIONS.md):

* fasta files whose basename ends in ``_notmerged.fas`` are dropped
  (shadow-pipeline artefacts — see ``[S04]``);
* empty fasta files are kept so empty samples reach the occurrence
  table (``[S09]``); downstream processes tolerate zero-record
  inputs;
* duplicate sample IDs abort the workflow (see ``[S13]`` /
  ``[S14]``).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import NamedTuple, Optional


# Suffix used by the Part A shadow pipeline ([S04]). Stripping the
# `.fas` extension first, anything ending in this token is dropped.
NOTMERGED_SUFFIX: str = "_notmerged"


class FastaSample(NamedTuple):
    """One Part B input fasta as seen by the workflow."""

    sample_id: str
    fasta: Path


class DuplicateSampleIDError(ValueError):
    """Raised when two or more fasta files derive the same sample ID."""


def _sample_id_from_fasta(path: Path) -> str:
    """Return the basename stripped of the ``.fas`` extension."""
    name = path.name
    if name.endswith(".fas"):
        return name[: -len(".fas")]
    return name


def _is_notmerged(path: Path) -> bool:
    return _sample_id_from_fasta(path).endswith(NOTMERGED_SUFFIX)


def discover(folders: list[Path]) -> list[FastaSample]:
    """Walk ``folders`` (no recursion) and return the Part B fasta channel.

    Files whose basename ends in ``_notmerged.fas`` are skipped.
    Empty files are kept so empty samples travel through to the
    occurrence table ([S09]).
    """
    samples: list[FastaSample] = []
    for folder in folders:
        for path in sorted(folder.glob("*.fas")):
            if not path.is_file():
                continue
            if _is_notmerged(path):
                continue
            samples.append(
                FastaSample(
                    sample_id=_sample_id_from_fasta(path), fasta=path
                )
            )
    return samples


def check_unique_sample_ids(samples: list[FastaSample]) -> None:
    """Reject any duplicate sample ID across the discovered fastas.

    Mirrors :func:`discover_fastq.check_unique_sample_ids`; the error
    message groups offending fasta paths per sample ID.
    """
    groups: dict[str, list[FastaSample]] = {}
    for sample in samples:
        groups.setdefault(sample.sample_id, []).append(sample)

    duplicates = {
        sid: members for sid, members in groups.items() if len(members) > 1
    }
    if not duplicates:
        return

    lines: list[str] = [
        "duplicate sample IDs (each sample ID must be unique):"
    ]
    for sid in sorted(duplicates):
        paths = [str(member.fasta) for member in duplicates[sid]]
        lines.append(f"  {sid}: {', '.join(paths)}")
    raise DuplicateSampleIDError("\n".join(lines))


# ---------- CLI ------------------------------------------------------------

def _parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Discover .fas files in one or more folders and emit a TSV "
            "(sample_id, fasta_path)."
        ),
    )
    parser.add_argument(
        "folders",
        nargs="+",
        type=Path,
        help="Folders to scan (no recursion).",
    )
    return parser.parse_args(argv)


def main(argv: Optional[list[str]] = None) -> int:
    args = _parse_args(argv)
    samples = discover(args.folders)
    try:
        check_unique_sample_ids(samples)
    except DuplicateSampleIDError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    for sample in samples:
        print(f"{sample.sample_id}\t{sample.fasta}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
