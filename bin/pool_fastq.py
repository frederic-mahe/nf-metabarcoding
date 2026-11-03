#!/usr/bin/env python3
"""Pool several fastq files into one gzipped fastq file ([S129]).

Files are concatenated in the order given (Part A sorts them by path
before calling this helper):

* ``.gz`` files are copied byte for byte — concatenated gzip members
  form a valid multi-member gzip file;
* plain and ``.bz2`` files are recompressed with gzip, with no file
  name and a timestamp of 0, so the output is byte-identical across
  runs.

The pooling log starts with ``<sample>: <N> files pooled`` and lists
one input file name per line. Used by the ``pool_fastq`` process for
folder pooling (``--pool_fastq_by_folder``, ``[S128]``) and for
samplesheet pooling (``[S132]``).
"""

import sys

if sys.version_info < (3, 10):
    sys.exit(
        "pool_fastq.py: Python >= 3.10 required, found "
        + sys.version.split()[0]
        + ". Enable a container or conda profile"
        + " (e.g. -profile apptainer)."
    )

import argparse  # noqa: E402
import bz2  # noqa: E402
import gzip  # noqa: E402
import shutil  # noqa: E402
from pathlib import Path  # noqa: E402
from typing import BinaryIO, Optional  # noqa: E402


# Copy buffer size: large enough to stream multi-GB inputs efficiently.
CHUNK: int = 1024 * 1024


def _recompress(source: BinaryIO, out: BinaryIO) -> None:
    """Append ``source``'s bytes to ``out`` as one gzip member."""
    with gzip.GzipFile(
        filename="", mode="wb", fileobj=out, mtime=0
    ) as member:
        shutil.copyfileobj(source, member, CHUNK)


def append(path: Path, out: BinaryIO) -> None:
    """Append one fastq file to the gzip stream ``out``."""
    if path.name.endswith(".gz"):
        with open(path, "rb") as source:
            shutil.copyfileobj(source, out, CHUNK)
    elif path.name.endswith(".bz2"):
        with bz2.open(path, "rb") as source:
            _recompress(source, out)
    else:
        with open(path, "rb") as source:
            _recompress(source, out)


def pool(sample: str, files: list[Path], output: Path, log: Path) -> None:
    """Concatenate ``files`` into ``output`` and write the pooling log."""
    with open(output, "wb") as out:
        for path in files:
            append(path, out)
    lines = [f"{sample}: {len(files)} files pooled"]
    lines.extend(path.name for path in files)
    log.write_text("\n".join(lines) + "\n")


def _parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Pool fastq files into one gzipped fastq file.",
    )
    parser.add_argument("--sample", required=True,
                        help="sample ID (first word of the log)")
    parser.add_argument("--output", required=True, type=Path,
                        help="pooled gzipped fastq file to write")
    parser.add_argument("--log", required=True, type=Path,
                        help="pooling log to write")
    parser.add_argument("files", nargs="+", type=Path,
                        help="fastq files to pool, in order")
    return parser.parse_args(argv)


def main(argv: Optional[list[str]] = None) -> int:
    args = _parse_args(argv)
    pool(args.sample, args.files, args.output, args.log)
    return 0


if __name__ == "__main__":
    sys.exit(main())
