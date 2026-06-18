#!/usr/bin/env python3
"""Validate and normalize the ``--input`` samplesheet ([S70]).

The samplesheet is a header-keyed CSV in one of two profiles, chosen by
the columns present:

* **fastq** (Part A): ``sample``, ``fastq_1``, ``fastq_2`` (optional),
  ``run`` (optional);
* **fasta** (Part B standalone): ``sample``, ``fasta``, ``qual``
  (optional), ``stats`` (optional).

This module performs **structural** validation only — columns, required
cells, duplicate / reserved sample IDs, single-end inference, and
sibling ``qual`` / ``stats`` defaulting. File existence is enforced
downstream by Nextflow's ``file(checkIfExists: true)`` when each row is
staged, so the helper needs no filesystem access and is fully
unit-testable from strings. Path cells get shell ``~`` expansion
([S60]); relative paths are left as-is for Nextflow to resolve against
the launch directory.

CLI: ``parse_samplesheet.py SHEET.csv`` prints the normalized rows as
TSV to stdout (the workflow consumes them via ``splitCsv``); structural
errors print to stderr and exit non-zero. The normalized header carries
the profile implicitly (``fastq_1`` vs ``fasta``), so the workflow
infers the profile the same way this module does.
"""

from __future__ import annotations

import csv
import os
import sys
from dataclasses import astuple, dataclass
from pathlib import Path
from typing import Union

RESERVED_SUFFIX = "notmerged"

FASTQ_REQUIRED = ("sample", "fastq_1")
FASTQ_OPTIONAL = ("fastq_2", "run")
FASTA_REQUIRED = ("sample", "fasta")
FASTA_OPTIONAL = ("qual", "stats")


class SamplesheetError(Exception):
    """Structural validation failure in the ``--input`` samplesheet."""


@dataclass
class FastqSample:
    sample: str
    fastq_1: str
    fastq_2: str  # empty string marks a single-end sample
    run: str

    @property
    def single_end(self) -> bool:
        return not self.fastq_2


@dataclass
class FastaSample:
    sample: str
    fasta: str
    qual: str
    stats: str

    @property
    def shadow(self) -> bool:
        return self.sample.endswith(RESERVED_SUFFIX)


Record = Union[FastqSample, FastaSample]


def _expand(path: str) -> str:
    """Apply shell ``~`` expansion ([S60]); leave the rest untouched."""
    return os.path.expanduser(path) if path else path


def infer_profile(fieldnames: list[str]) -> str:
    """Return ``'fastq'`` or ``'fasta'`` from the header columns."""
    if "fastq_1" in fieldnames:
        return "fastq"
    if "fasta" in fieldnames:
        return "fasta"
    raise SamplesheetError(
        "cannot infer samplesheet profile: header has neither a "
        "'fastq_1' column (fastq profile) nor a 'fasta' column "
        f"(fasta profile); got columns {fieldnames}"
    )


def _check_columns(
    fieldnames: list[str], required: tuple, optional: tuple
) -> None:
    allowed = set(required) | set(optional)
    missing = [c for c in required if c not in fieldnames]
    if missing:
        raise SamplesheetError(
            f"samplesheet header is missing required column(s): {missing}"
        )
    unknown = [c for c in fieldnames if c not in allowed]
    if unknown:
        raise SamplesheetError(
            f"samplesheet header has unknown column(s): {unknown}; "
            f"allowed columns for this profile are {sorted(allowed)}"
        )


def _require(row: dict[str, str], column: str, line: int) -> str:
    value = row.get(column, "").strip()
    if not value:
        raise SamplesheetError(
            f"row {line}: required column '{column}' is empty"
        )
    return value


def _check_unique(records: list[Record], lines: list[int]) -> None:
    seen: dict[str, list[int]] = {}
    for record, line in zip(records, lines):
        seen.setdefault(record.sample, []).append(line)
    dupes = {s: ls for s, ls in seen.items() if len(ls) > 1}
    if dupes:
        groups = "; ".join(
            f"{sample} (rows {', '.join(map(str, ls))})"
            for sample, ls in dupes.items()
        )
        raise SamplesheetError(
            f"duplicate sample IDs (each must be unique): {groups}"
        )


def _parse_fastq(rows: list[tuple[int, dict]]) -> list[FastqSample]:
    records: list[FastqSample] = []
    lines: list[int] = []
    for line, row in rows:
        sample = _require(row, "sample", line)
        if sample.endswith(RESERVED_SUFFIX):
            raise SamplesheetError(
                f"row {line}: sample '{sample}' ends with the reserved "
                f"suffix '{RESERVED_SUFFIX}' (column 'sample')"
            )
        records.append(
            FastqSample(
                sample=sample,
                fastq_1=_expand(_require(row, "fastq_1", line)),
                fastq_2=_expand(row.get("fastq_2", "").strip()),
                run=row.get("run", "").strip(),
            )
        )
        lines.append(line)
    _check_unique(records, lines)
    return records


def _parse_fasta(rows: list[tuple[int, dict]]) -> list[FastaSample]:
    records: list[FastaSample] = []
    lines: list[int] = []
    for line, row in rows:
        sample = _require(row, "sample", line)
        fasta = _expand(_require(row, "fasta", line))
        qual = row.get("qual", "").strip()
        stats = row.get("stats", "").strip()
        parent = Path(fasta).parent
        records.append(
            FastaSample(
                sample=sample,
                fasta=fasta,
                qual=_expand(qual) if qual else str(parent / f"{sample}.qual"),
                stats=(
                    _expand(stats)
                    if stats
                    else str(parent / f"{sample}.stats")
                ),
            )
        )
        lines.append(line)
    _check_unique(records, lines)
    return records


def parse_rows(
    fieldnames: list[str], rows: list[tuple[int, dict]]
) -> tuple[str, list[Record]]:
    """Validate header + rows; return ``(profile, records)``."""
    profile = infer_profile(fieldnames)
    if profile == "fastq":
        _check_columns(fieldnames, FASTQ_REQUIRED, FASTQ_OPTIONAL)
        return profile, _parse_fastq(rows)
    _check_columns(fieldnames, FASTA_REQUIRED, FASTA_OPTIONAL)
    return profile, _parse_fasta(rows)


def read_samplesheet(
    path: Union[str, Path]
) -> tuple[list[str], list[tuple[int, dict]]]:
    """Read a CSV into ``(fieldnames, [(line_number, row_dict), ...])``."""
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle)
        try:
            header = next(reader)
        except StopIteration:
            raise SamplesheetError("samplesheet is empty (no header row)")
        fieldnames = [cell.strip() for cell in header]
        rows: list[tuple[int, dict]] = []
        for line, raw in enumerate(reader, start=2):
            if not any(cell.strip() for cell in raw):
                continue  # skip blank lines
            cells = [cell.strip() for cell in raw]
            cells += [""] * (len(fieldnames) - len(cells))
            rows.append((line, dict(zip(fieldnames, cells))))
    return fieldnames, rows


def to_tsv(profile: str, records: list[Record]) -> str:
    """Render normalized records as TSV (header + one row each)."""
    if profile == "fastq":
        header = ("sample", "fastq_1", "fastq_2", "run")
    else:
        header = ("sample", "fasta", "qual", "stats")
    out = ["\t".join(header)]
    for record in records:
        out.append("\t".join(astuple(record)))
    return "\n".join(out) + "\n"


def main(argv: Union[list[str], None] = None) -> int:
    args = sys.argv[1:] if argv is None else argv
    if len(args) != 1:
        sys.stderr.write("usage: parse_samplesheet.py SHEET.csv\n")
        return 2
    try:
        fieldnames, rows = read_samplesheet(args[0])
        profile, records = parse_rows(fieldnames, rows)
    except SamplesheetError as err:
        sys.stderr.write(f"error: {err}\n")
        return 1
    sys.stdout.write(to_tsv(profile, records))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
