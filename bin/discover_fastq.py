#!/usr/bin/env python3
"""Discover fastq files in one or more folders and identify R1/R2 pairs.

Single source of truth for the canonical paired-end name patterns
documented in ``SPECIFICATIONS.md`` (the "Common fastq file-name
patterns" table). Both R2-name derivation and sample-ID extraction
walk the same ``PATTERNS`` list, so adding a new row updates both
operations at once.

Used as:

* an importable module — call :func:`derive_r2_name`,
  :func:`derive_sample_id`, or :func:`discover` directly;
* a CLI — ``discover_fastq.py FOLDER [FOLDER ...]`` emits TSV
  ``sample_id\\tr1\\tr2`` (R2 is empty for single-end samples) to
  stdout.

See ``[S10]``, ``[S11]``, ``[S12]``, ``[S21]`` in SPECIFICATIONS.md.
"""

from __future__ import annotations

import argparse
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import NamedTuple, Optional


# ---------- pattern table --------------------------------------------------

# Every regex anchors to the end of the basename via ``$``. Named
# groups: ``sample`` (prefix that becomes the sample ID) and any
# others referenced by ``r2_template``. ``r2_template`` is a
# ``str.format``-style string whose placeholders are the named groups
# of ``regex`` — this is how R1→R2 reuses the same pattern definition.
@dataclass(frozen=True)
class PatternEntry:
    """One canonical naming convention.

    ``regex`` matches an R1 basename. ``r2_template`` reconstructs the
    paired R2 basename from the named groups in ``regex``.
    """

    name: str
    regex: re.Pattern[str]
    r2_template: str


# Shared sub-expression for fastq extensions. The named group is
# referenced by ``r2_template`` so every pattern can reuse it without
# duplicating the alternation.
_EXT = r"(?P<ext>\.(?:fastq|fq)(?:\.(?:gz|bz2))?)"

PATTERNS: list[PatternEntry] = [
    # row 1 — MiSeq default: `_L00[1-9]_R1_00[1-9].<ext>`
    PatternEntry(
        name="miseq_001",
        regex=re.compile(
            rf"^(?P<sample>.*)(?P<lane>_L00[1-9])_R1(?P<tail>_00[1-9]){_EXT}$"
        ),
        r2_template="{sample}{lane}_R2{tail}{ext}",
    ),
    # row 2 — MiSeq with middle segment: `_L00[1-9]_<mid>_R1.<ext>`
    PatternEntry(
        name="miseq_mid",
        regex=re.compile(
            rf"^(?P<sample>.*)(?P<lane>_L00[1-9])_(?P<mid>.+)_R1{_EXT}$"
        ),
        r2_template="{sample}{lane}_{mid}_R2{ext}",
    ),
    # row 3 — plain MiSeq: `_L00[1-9]_R1.<ext>`
    PatternEntry(
        name="miseq_plain",
        regex=re.compile(
            rf"^(?P<sample>.*)(?P<lane>_L00[1-9])_R1{_EXT}$"
        ),
        r2_template="{sample}{lane}_R2{ext}",
    ),
    # row 4 — numeric-lane with tail: `[._][1-9]_1_<tail>.<ext>`
    PatternEntry(
        name="numeric_lane_tail",
        regex=re.compile(
            rf"^(?P<sample>.*)(?P<sep>[._])(?P<lane>[1-9])_1_"
            rf"(?P<tail>[^/]+?){_EXT}$"
        ),
        r2_template="{sample}{sep}{lane}_2_{tail}{ext}",
    ),
    # row 5 — numeric-lane: `[._][1-9]_1.<ext>`
    PatternEntry(
        name="numeric_lane",
        regex=re.compile(
            rf"^(?P<sample>.*)(?P<sep>[._])(?P<lane>[1-9])_1{_EXT}$"
        ),
        r2_template="{sample}{sep}{lane}_2{ext}",
    ),
    # row 6 — generic R1: `[._]R1.<ext>`
    PatternEntry(
        name="generic_r1",
        regex=re.compile(
            rf"^(?P<sample>.*)(?P<sep>[._])R1{_EXT}$"
        ),
        r2_template="{sample}{sep}R2{ext}",
    ),
    # row 7 — generic 1: `[._]1.<ext>`
    PatternEntry(
        name="generic_1",
        regex=re.compile(
            rf"^(?P<sample>.*)(?P<sep>[._])1{_EXT}$"
        ),
        r2_template="{sample}{sep}2{ext}",
    ),
]


# Acceptable extensions for the discovery glob; mirrors `_EXT` above.
_FASTQ_GLOBS: tuple[str, ...] = (
    "*.fastq",
    "*.fq",
    "*.fastq.gz",
    "*.fq.gz",
    "*.fastq.bz2",
    "*.fq.bz2",
)


# ---------- user --fastq_pattern → PatternEntry ----------------------------

# Translates the glob meta-chars users actually need. We deliberately
# keep the surface small: `*` (any-chars), `{R1,R2}` / `{1,2}` (the
# discriminator), literal `.`. Anything else is treated literally,
# which matches the meaning the workflow had before this change.
def _user_pattern_to_entry(glob: str) -> PatternEntry:
    """Convert a user-supplied glob into a :class:`PatternEntry`.

    The glob must contain exactly one ``{<r1>,<r2>}`` brace token
    (commonly ``{1,2}`` or ``{R1,R2}``); this is the R1/R2
    discriminator. The first ``*`` before the discriminator becomes
    the sample-ID capture group.
    """
    brace = re.search(r"\{([^,{}]+),([^,{}]+)\}", glob)
    if brace is None:
        raise ValueError(
            f"--fastq_pattern must contain a brace token like {{1,2}} "
            f"or {{R1,R2}} marking the R1/R2 discriminator; got {glob!r}"
        )

    r1_token, r2_token = brace.group(1), brace.group(2)
    prefix, suffix = glob[: brace.start()], glob[brace.end():]

    sample_taken = False

    def _translate(segment: str, *, capture_sample: bool) -> str:
        nonlocal sample_taken
        out: list[str] = []
        for ch in segment:
            if ch == "*":
                if capture_sample and not sample_taken:
                    out.append("(?P<sample>.*)")
                    sample_taken = True
                else:
                    out.append(".*")
            else:
                out.append(re.escape(ch))
        return "".join(out)

    prefix_re = _translate(prefix, capture_sample=True)
    suffix_re = _translate(suffix, capture_sample=False)

    if not sample_taken:
        # No `*` in the prefix — anchor sample to the empty match so
        # the named group exists and ``derive_sample_id`` still works.
        prefix_re = "(?P<sample>)" + prefix_re

    r1_regex = re.compile(f"^{prefix_re}{re.escape(r1_token)}{suffix_re}$")
    # Re-quote braces so str.format leaves them alone in the literal text.
    template_prefix = prefix.replace("{", "{{").replace("}", "}}")
    template_suffix = suffix.replace("{", "{{").replace("}", "}}")
    # The sample-ID `*` placeholder becomes a {sample} format slot; any
    # other `*` characters become `*` literals (rare in practice).
    template_prefix_filled = template_prefix.replace("*", "{sample}", 1)
    r2_template = f"{template_prefix_filled}{r2_token}{template_suffix}"

    return PatternEntry(name="user", regex=r1_regex, r2_template=r2_template)


# ---------- public API -----------------------------------------------------

class Sample(NamedTuple):
    """One sample as seen by the workflow."""

    sample_id: str
    r1: Path
    r2: Optional[Path]  # None for single-end / unpaired samples


def _walk(
    name: str, entries: list[PatternEntry]
) -> Optional[tuple[PatternEntry, re.Match[str]]]:
    """Return the first ``(entry, match)`` pair that matches ``name``."""
    for entry in entries:
        m = entry.regex.match(name)
        if m is not None:
            return entry, m
    return None


def derive_r2_name(
    r1_name: str, *, extra_pattern: Optional[str] = None
) -> Optional[str]:
    """Return the R2 basename paired with ``r1_name``, or ``None``.

    ``extra_pattern`` (a user-supplied glob) is checked **before** the
    canonical patterns, mirroring the precedence rule in
    ``[S11]``.
    """
    entries = _entries_with_optional_override(extra_pattern)
    hit = _walk(r1_name, entries)
    if hit is None:
        return None
    entry, m = hit
    return entry.r2_template.format(**m.groupdict())


def derive_sample_id(
    r1_name: str, *, extra_pattern: Optional[str] = None
) -> Optional[str]:
    """Return the sample ID for ``r1_name``, or ``None``.

    The sample ID is the ``sample`` capture group of the first
    matching pattern (``[S12]``).
    """
    entries = _entries_with_optional_override(extra_pattern)
    hit = _walk(r1_name, entries)
    if hit is None:
        return None
    return hit[1]["sample"]


def _entries_with_optional_override(
    extra_pattern: Optional[str],
) -> list[PatternEntry]:
    if extra_pattern is None:
        return PATTERNS
    return [_user_pattern_to_entry(extra_pattern), *PATTERNS]


def _strip_extension(name: str) -> str:
    """Strip the ``.(fastq|fq)(.(gz|bz2))?`` suffix from ``name``."""
    return re.sub(r"\.(fastq|fq)(\.(gz|bz2))?$", "", name)


def _collect_fastq(folders: list[Path]) -> list[Path]:
    """Glob every fastq file in the listed folders (no recursion)."""
    seen: list[Path] = []
    for folder in folders:
        for pattern in _FASTQ_GLOBS:
            seen.extend(sorted(folder.glob(pattern)))
    return seen


def discover(
    folders: list[Path], *, extra_pattern: Optional[str] = None
) -> list[Sample]:
    """Walk ``folders`` and return the list of discovered samples.

    Pairing rules:
    * an R1 whose R2 partner exists in the same folder becomes a
      paired-end ``Sample`` (sample ID per ``[S12]``);
    * an R1 whose R2 partner is missing falls through to single-end
      (decision-4a, "forgiving");
    * any file that matches no R1 pattern is single-end (``[S21]``);
    * an R2 file consumed by its R1 partner is NOT re-emitted.
    """
    all_paths = _collect_fastq(folders)
    by_path: dict[Path, Path] = {p: p for p in all_paths}
    consumed_as_r2: set[Path] = set()
    samples: list[Sample] = []

    for path in all_paths:
        if path in consumed_as_r2:
            continue

        sample_id = derive_sample_id(path.name, extra_pattern=extra_pattern)
        r2_name = derive_r2_name(path.name, extra_pattern=extra_pattern)

        if sample_id is not None and r2_name is not None:
            r2_path = path.with_name(r2_name)
            if r2_path in by_path:
                samples.append(
                    Sample(sample_id=sample_id, r1=path, r2=r2_path)
                )
                consumed_as_r2.add(r2_path)
                continue
            # Orphan R1 — fall through to single-end.
            print(
                f"warning: {path.name} matches a paired-end pattern but "
                f"its R2 partner {r2_name} is missing; treating as single-end",
                file=sys.stderr,
            )

        samples.append(
            Sample(sample_id=_strip_extension(path.name), r1=path, r2=None)
        )

    return samples


# ---------- CLI ------------------------------------------------------------

def _parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Discover fastq files in one or more folders and emit a TSV "
            "(sample_id, r1, r2) — r2 is empty for single-end samples."
        ),
    )
    parser.add_argument(
        "folders",
        nargs="+",
        type=Path,
        help="Folders to scan (no recursion).",
    )
    parser.add_argument(
        "--extra-pattern",
        type=str,
        default=None,
        help=(
            "User-supplied paired-end glob with a brace token marking the "
            "R1/R2 discriminator (e.g. *_run17_{1,2}.fastq.gz). Checked "
            "before the canonical pattern table."
        ),
    )
    return parser.parse_args(argv)


def main(argv: Optional[list[str]] = None) -> int:
    args = _parse_args(argv)
    for sample in discover(args.folders, extra_pattern=args.extra_pattern):
        r2 = str(sample.r2) if sample.r2 is not None else ""
        print(f"{sample.sample_id}\t{sample.r1}\t{r2}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
