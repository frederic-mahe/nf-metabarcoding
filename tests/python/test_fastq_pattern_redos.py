"""Cap the wildcard count in ``--fastq_pattern`` ([S96]).

Each ``*`` in a user pattern becomes a greedy ``.*`` in the discovery
regex. An unbounded number of them lets a crafted pattern drive
catastrophic regex backtracking — a denial-of-service hang — against the
scanned file names. A legitimate pattern needs one ``*`` to capture the
sample prefix (occasionally a second), so a small fixed cap never blocks
a real pattern while removing the pathological case. The cap fires at
pattern-parse time, before any file is globbed.
"""

# COVERAGE: [S96]

from __future__ import annotations

import pytest

from discover_fastq import _user_pattern_to_entry


# ---------- patterns within the cap are accepted --------------------------

@pytest.mark.parametrize(
    "glob",
    [
        "*_{1,2}.fastq.gz",          # one wildcard (the common case)
        "*_*_{1,2}.fastq.gz",        # two wildcards
        "*_*_*_{R1,R2}.fastq.gz",    # three wildcards (the cap)
    ],
)
def test_wildcard_count_within_cap_accepted(glob: str) -> None:
    # no exception
    _user_pattern_to_entry(glob)


# ---------- patterns over the cap are rejected ----------------------------

@pytest.mark.parametrize(
    "glob",
    [
        "*a*a*a*a{1,2}.fastq",           # four wildcards
        "*_*_*_*_*_{1,2}.fastq.gz",      # five wildcards
    ],
)
def test_wildcard_count_over_cap_rejected(glob: str) -> None:
    with pytest.raises(ValueError) as excinfo:
        _user_pattern_to_entry(glob)
    # the message states the limit so the user can fix the pattern
    assert "3" in str(excinfo.value)
