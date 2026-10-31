"""Part C taxonomy helpers pass non-UTF-8 header bytes through unchanged."""

# COVERAGE: [S112]

import gzip
import os
import subprocess
import sys
from pathlib import Path

import pytest

BIN = Path(__file__).resolve().parents[2] / "bin"

# a Latin-1 'é' (not valid UTF-8) and a UTF-8 'é'
LINEAGE = b"Eukaryota|Genus_\xe9|Species_\xc3\xa9"

HEADER = (
    b"OTU\ttotal\tcloud\tamplicon\tlength\tabundance\tchimera\tspread\t"
    b"quality\tsequence\tidentity\ttaxonomy\treferences\ts1\n"
)


def _run(script: str, *args: str) -> bytes:
    # LC_ALL=C: the bare-cluster-node case; the helpers must not depend
    # on the locale's encoding.
    env = {**os.environ, "LC_ALL": "C"}
    env.pop("PYTHONUTF8", None)
    result = subprocess.run(
        [sys.executable, str(BIN / script), *args],
        capture_output=True,
        env=env,
    )
    assert result.returncode == 0, result.stderr.decode(errors="replace")
    return result.stdout


def test_stampa_merge_passes_lineage_bytes(tmp_path: Path) -> None:
    hits = tmp_path / "hits"
    hits.write_bytes(b"amp1;size=10;\t100.0\trefA " + LINEAGE + b"\n")
    out = _run("stampa_merge.py", str(hits))
    assert out == b"amp1\t10\t100.0\t" + LINEAGE + b"\trefA\n"


def test_update_occurrence_table_passes_taxonomy_bytes(
    tmp_path: Path,
) -> None:
    table = tmp_path / "table.tsv"
    table.write_bytes(
        HEADER
        + b"1\t10\tNA\tamp1\t8\t10\tN\t1\t0.001\tACGTACGT\t0.0\tNA\tNA\t10\n"
    )
    assignments = tmp_path / "assignments.tsv"
    assignments.write_bytes(b"amp1\t10\t100.0\t" + LINEAGE + b"\trefA\n")
    out = _run(
        "update_occurrence_table.py",
        "--occurrence_table", str(table),
        "--assignments", str(assignments),
    )
    row = out.splitlines()[1].split(b"\t")
    assert row[11] == LINEAGE


@pytest.mark.parametrize("compressed", [False, True])
def test_majority_assignment_passes_reference_bytes(
    tmp_path: Path, compressed: bool
) -> None:
    fasta = b">refA " + LINEAGE + b"\nACGTACGT\n"
    if compressed:
        reference = tmp_path / "ref.fas.gz"
        reference.write_bytes(gzip.compress(fasta))
    else:
        reference = tmp_path / "ref.fas"
        reference.write_bytes(fasta)
    table = tmp_path / "table_assigned.tsv"
    table.write_bytes(
        HEADER
        + b"1\t10\tNA\tamp1\t8\t10\tN\t1\t0.001\tACGTACGT\t100.0\t"
        + LINEAGE + b"\trefA\t10\n"
    )
    out = _run(
        "majority_assignment.py",
        "--input_table", str(table),
        "--reference_db", str(reference),
    )
    assert out.splitlines()[1] == (
        b"1\tamp1\tEukaryota (1/1)|Genus_\xe9 (1/1)|Species_\xc3\xa9 (1/1)"
    )
