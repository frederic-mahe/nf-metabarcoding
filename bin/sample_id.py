#!/usr/bin/env python3
"""Shared sample-ID charset validator ([S93]).

A sample ID flows verbatim into the Part A / Part B Nextflow process
scripts as a shell token and an output-file basename (``!{sampleId}``,
e.g. ``--log !{sampleId}_merging.log``). Nextflow substitutes the raw
string into ``.command.sh`` before bash parses it, so an ID carrying a
shell metacharacter (``$(...)``, backticks, ``;``, ``|`` …) would be
executed, and one carrying whitespace or a path separator would break
or escape the output filename.

To close that off at the two entry points at once — the ``--input``
samplesheet (:mod:`parse_samplesheet`) and folder discovery
(:mod:`discover_fastq` / :mod:`discover_fasta`) — both import
:func:`validate_sample_id` from here. Keeping the rule in one module
means the samplesheet and discovery paths cannot drift apart.

The safe set is ``[A-Za-z0-9._-]`` with the extra rule that the first
character must be a letter, digit, or underscore: a leading ``-`` would
be read as an option by downstream tools and a leading ``.`` would name
a hidden file or open ``..`` traversal.
"""

from __future__ import annotations

import re

# First character: letter, digit, or underscore (no leading `-` or `.`).
# Remaining characters: the safe set plus `.` and `-`.
_VALID_SAMPLE_ID = re.compile(r"^[A-Za-z0-9_][A-Za-z0-9._-]*$")


class InvalidSampleIdError(ValueError):
    """Raised when a sample ID falls outside the safe character set."""


def is_valid_sample_id(sample_id: str) -> bool:
    """Return ``True`` when ``sample_id`` is safe to use as a shell token."""
    return bool(_VALID_SAMPLE_ID.match(sample_id))


def validate_sample_id(sample_id: str, *, context: str = "") -> str:
    """Return ``sample_id`` unchanged when valid; raise otherwise.

    ``context`` (e.g. ``"row 4"``) is prepended to the error message so
    the caller can point at the offending samplesheet row or file.
    """
    if is_valid_sample_id(sample_id):
        return sample_id
    where = f"{context}: " if context else ""
    raise InvalidSampleIdError(
        f"{where}sample ID {sample_id!r} is not allowed; a sample ID must "
        "start with a letter, digit, or underscore and contain only "
        "letters, digits, '.', '_', and '-' (no whitespace, path "
        "separators, or shell metacharacters)"
    )
