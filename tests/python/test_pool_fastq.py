"""Tests for ``bin/pool_fastq.py`` ([S129]).

The helper concatenates several fastq files, in the order given, into
one gzip file: ``.gz`` inputs are copied byte for byte, plain and
``.bz2`` inputs are recompressed with gzip (no name, timestamp 0), so
the output is byte-identical across runs. It also writes a pooling log.
"""

# COVERAGE: [S129]

import bz2
import gzip
import subprocess
import sys
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[2] / "bin" / "pool_fastq.py"

RECORDS = [
    b"@r0\nACGT\n+\nIIII\n",
    b"@r1\nCCGG\n+\nIIII\n",
    b"@r2\nTTAA\n+\nIIII\n",
]


def _inputs(folder: Path) -> list[Path]:
    folder.mkdir(parents=True, exist_ok=True)
    gz = folder / "a_0.fastq.gz"
    gz.write_bytes(gzip.compress(RECORDS[0], mtime=0))
    plain = folder / "b_1.fastq"
    plain.write_bytes(RECORDS[1])
    bz = folder / "c_2.fq.bz2"
    bz.write_bytes(bz2.compress(RECORDS[2]))
    return [gz, plain, bz]


def _run(sample: str, out: Path, log: Path, files: list[Path]
         ) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, str(SCRIPT),
         "--sample", sample,
         "--output", str(out),
         "--log", str(log),
         *[str(f) for f in files]],
        capture_output=True,
        text=True,
    )


def test_pools_gz_plain_and_bz2_in_given_order(tmp_path: Path) -> None:
    files = _inputs(tmp_path / "in")
    out = tmp_path / "S.fastq.gz"
    log = tmp_path / "S_pooling.log"

    proc = _run("S", out, log, files)

    assert proc.returncode == 0, proc.stderr
    assert gzip.decompress(out.read_bytes()) == b"".join(RECORDS)

    reversed_out = tmp_path / "R.fastq.gz"
    proc = _run("R", reversed_out, tmp_path / "R.log", files[::-1])
    assert proc.returncode == 0, proc.stderr
    assert gzip.decompress(reversed_out.read_bytes()) == (
        b"".join(RECORDS[::-1])
    )


def test_gz_members_are_copied_byte_for_byte(tmp_path: Path) -> None:
    files = _inputs(tmp_path / "in")
    out = tmp_path / "S.fastq.gz"

    _run("S", out, tmp_path / "S.log", files)

    assert out.read_bytes().startswith(files[0].read_bytes())


def test_output_is_byte_identical_across_runs(tmp_path: Path) -> None:
    files = _inputs(tmp_path / "in")
    first = tmp_path / "first.fastq.gz"
    second = tmp_path / "second.fastq.gz"

    _run("S", first, tmp_path / "1.log", files)
    _run("S", second, tmp_path / "2.log", files)

    assert first.read_bytes() == second.read_bytes()


def test_log_counts_and_lists_files(tmp_path: Path) -> None:
    files = _inputs(tmp_path / "in")
    log = tmp_path / "S_pooling.log"

    _run("S", tmp_path / "S.fastq.gz", log, files)

    lines = log.read_text().splitlines()
    assert lines[0] == "S: 3 files pooled", lines
    assert lines[1:] == [f.name for f in files], lines


def test_empty_inputs_give_an_empty_fastq(tmp_path: Path) -> None:
    folder = tmp_path / "in"
    folder.mkdir()
    a = folder / "a_0.fastq"
    a.write_bytes(b"")
    b = folder / "b_1.fastq.gz"
    b.write_bytes(gzip.compress(b"", mtime=0))
    out = tmp_path / "S.fastq.gz"

    proc = _run("S", out, tmp_path / "S.log", [a, b])

    assert proc.returncode == 0, proc.stderr
    assert gzip.decompress(out.read_bytes()) == b""
