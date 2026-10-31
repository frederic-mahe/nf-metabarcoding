"""Empty ranks in a stampa reference lineage (`A||C|D`)."""

# COVERAGE: [S113]

import io
from pathlib import Path

import pytest

from majority_assignment import majority_for_references
from stampa_merge import emit_results


@pytest.mark.parametrize(
    ("lineages", "expected"),
    [
        # every top hit carries the empty rank: kept as is
        (["A||C|D", "A||C|D"], "A||C|D"),
        # empty vs named rank: a disagreement like any other
        (["A||C|D", "A|B|C|D"], "A|*|C|D"),
        # a trailing empty rank is kept too
        (["A|B|C|", "A|B|C|"], "A|B|C|"),
    ],
)
def test_stampa_merge_keeps_empty_ranks(
    tmp_path: Path, lineages: list[str], expected: str
) -> None:
    hits = tmp_path / "hits"
    hits.write_text(
        "".join(
            f"q1;size=10\t100.0\tr{i} {lineage}\n"
            for i, lineage in enumerate(lineages, start=1)
        )
    )
    out = io.StringIO()
    emit_results(str(hits), out)
    references = ",".join(f"r{i}" for i in range(1, len(lineages) + 1))
    assert out.getvalue().splitlines() == [
        f"q1\t10\t100.0\t{expected}\t{references}"
    ]


def test_majority_counts_empty_rank_as_a_name() -> None:
    # The empty name wins rank 2 (2 votes vs 1 for B) and is written
    # with no name; it is not replaced with NA.
    accessions = {"r1": "A||C|D", "r2": "A||C|D", "r3": "A|B|C|D"}
    assert majority_for_references("r1,r2,r3", accessions) == \
        "A (3/3)| (2/3)|C (3/3)|D (3/3)"
