#!/usr/bin/env python3
"""Check every header of a taxonomic reference before Part C assignment.

`[S114]`-`[S117]` (D22): the startup sniff (`[S73]`) reads only the first
header; this script reads them all, once, as an ordinary task. Each
header is checked against the grammar of the method that reads the
file:

- stampa: ``>accession<SPACE>lineage``, the accession being the text
  before the first whitespace and the lineage the ``|``-separated text
  after it;
- sintax: ``>accession;tax=rank:name,rank:name,...`` terminated by
  ``;`` or by the end of the header.

Findings are written to a TSV report (``line``, ``accession``,
``issue``) and summarised in a text file whose first line is the
status (``clean``, ``warning`` or ``error``). The script exits 0 even
when it finds errors: the workflow publishes the report, then stops the
run (`[S117]`).
"""

import sys

if sys.version_info < (3, 10):
    sys.exit(
        "check_reference.py: Python >= 3.10 required, found "
        + sys.version.split()[0]
        + ". Enable a container or conda profile"
        + " (e.g. -profile apptainer)."
    )

__author__ = "Frédéric Mahé <frederic.mahe@cirad.fr>"
__date__ = "2026/10/31"
__version__ = "$Revision: 1.0"

import argparse  # noqa: E402
import bz2  # noqa: E402
import gzip  # noqa: E402
import re  # noqa: E402
from collections import Counter  # noqa: E402
from typing import IO, Iterator, NamedTuple, Optional  # noqa: E402


# [S115]: abort the run before any assignment.
ERRORS = frozenset({
    "format",
    "invalid_utf8",
    "control_char",
    "duplicate_accession",
    "comma_in_accession",
    "whitespace",
})
# [S116]: reported, the run continues.
WARNINGS = frozenset({"empty_rank", "star_in_lineage"})

# gzip and bzip2 file signatures (first two bytes).
GZIP_MAGIC = b"\x1f\x8b"
BZIP2_MAGIC = b"BZ"

# number of offending headers quoted in the summary ([S117])
SUMMARY_HEADERS = 5

ACCESSION_RE = re.compile(r"[^;\s]*")
TAX_RE = re.compile(r";tax=([^;]*)")


class Finding(NamedTuple):
    line: int
    accession: str
    issue: str


def parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Check every header of a stampa or sintax reference and "
            "write a TSV report plus a status summary."
        ),
    )
    parser.add_argument(
        "--grammar", required=True, choices=("stampa", "sintax"),
        help="header grammar of the method that reads the reference",
    )
    parser.add_argument(
        "--parameter", required=True,
        help="parameter name quoted in the summary (e.g. reference_dataset)",
    )
    parser.add_argument(
        "--majority", action="store_true",
        help="reject commas in accessions (--majority_assignment)",
    )
    parser.add_argument("--report", required=True, help="TSV report path")
    parser.add_argument("--summary", required=True, help="summary path")
    parser.add_argument("reference", help="reference fasta (plain, .gz, .bz2)")
    return parser.parse_args(argv)


def open_binary(filename: str) -> IO[bytes]:
    """Open a reference as bytes, decompressing by magic bytes."""
    with open(filename, "rb") as probe:
        signature = probe.read(2)
    if signature == GZIP_MAGIC:
        return gzip.open(filename, "rb")
    if signature == BZIP2_MAGIC:
        return bz2.open(filename, "rb")
    return open(filename, "rb")


def iter_headers(filename: str) -> Iterator[tuple[int, bytes]]:
    """Yield ``(line number, header bytes)`` without ``>`` or line end.

    A carriage return ending the line (CRLF file) is dropped.
    """
    with open_binary(filename) as handle:
        for number, raw in enumerate(handle, start=1):
            if not raw.startswith(b">"):
                continue
            header = raw[1:].removesuffix(b"\n").removesuffix(b"\r")
            yield number, header


def has_control_char(text: str) -> bool:
    return any(ord(c) < 32 or ord(c) == 127 for c in text)


def stampa_issues(text: str) -> list[str]:
    """Grammar issues of a stampa header (accession, lineage)."""
    parts = re.split(r"\s", text, maxsplit=1)
    lineage = parts[1] if len(parts) == 2 else ""
    if not parts[0] or not lineage.strip():
        return ["format"]
    issues = []
    if lineage != lineage.strip():
        issues.append("whitespace")
    if "" in lineage.split("|"):
        issues.append("empty_rank")
    if "*" in lineage:
        issues.append("star_in_lineage")
    return issues


def sintax_issues(text: str, accession: str) -> list[str]:
    """Grammar issues of a sintax header (accession;tax=...)."""
    match = TAX_RE.search(text)
    if (
        not accession
        or match is None
        or not text.startswith(accession + ";")
        or not match.group(1)
    ):
        return ["format"]
    issues = []
    for field in match.group(1).split(","):
        if field == "":
            issues.append("empty_rank")
            continue
        rank, colon, name = field.partition(":")
        if not colon or not rank:
            return ["format"]
        if name == "":
            issues.append("empty_rank")
    # one finding per issue code and header
    return list(dict.fromkeys(issues))


def check_reference(
    filename: str, grammar: str, majority: bool = False
) -> list[Finding]:
    """Return every finding of `filename` under `grammar`, by line."""
    findings: list[Finding] = []
    seen: set[str] = set()
    for number, header in iter_headers(filename):
        issues = []
        try:
            text = header.decode("utf-8")
        except UnicodeDecodeError:
            issues.append("invalid_utf8")
            # [S112]: keep the bytes, so the report shows them as is
            text = header.decode("utf-8", errors="surrogateescape")
        if has_control_char(text):
            issues.append("control_char")
        if grammar == "stampa":
            accession = re.split(r"\s", text, maxsplit=1)[0]
            issues.extend(stampa_issues(text))
            if majority and "," in accession:
                issues.append("comma_in_accession")
        else:
            accession = ACCESSION_RE.match(text).group(0)
            issues.extend(sintax_issues(text, accession))
        if accession:
            if accession in seen:
                issues.append("duplicate_accession")
            seen.add(accession)
        findings.extend(Finding(number, accession, i) for i in issues)
    return findings


def status(findings: list[Finding]) -> str:
    """``error``, ``warning`` or ``clean``."""
    issues = {f.issue for f in findings}
    if issues & ERRORS:
        return "error"
    if issues:
        return "warning"
    return "clean"


def summarise(findings: list[Finding], parameter: str, report: str) -> str:
    """Counts per issue code, first offending headers, report path.

    Empty for a clean reference.
    """
    state = status(findings)
    if state == "clean":
        return ""
    if state == "error":
        intro = (
            f"--{parameter}: the reference has errors; the run stops "
            "before the taxonomic assignment."
        )
    else:
        intro = (
            f"--{parameter}: the reference has warnings; the taxonomic "
            "assignment runs anyway."
        )
    counts = Counter(f.issue for f in findings)
    lines = [intro, "Findings per issue:"]
    lines.extend(f"  {issue}: {n}" for issue, n in sorted(counts.items()))
    lines.append("First offending headers:")
    seen_lines: list[int] = []
    for finding in findings:
        if finding.line in seen_lines:
            continue
        if len(seen_lines) == SUMMARY_HEADERS:
            break
        seen_lines.append(finding.line)
        issues = ", ".join(
            f.issue for f in findings if f.line == finding.line
        )
        lines.append(
            f"  line {finding.line}: {finding.accession} ({issues})"
        )
    lines.append(f"Full report: {report}")
    return "\n".join(lines)


def write_report(findings: list[Finding], report: str) -> None:
    """Write the TSV report, header row first ([S117])."""
    # [S112]: write accession bytes back unchanged
    with open(
        report, "w", encoding="utf-8", errors="surrogateescape"
    ) as handle:
        print("line", "accession", "issue", sep="\t", file=handle)
        for finding in findings:
            print(*finding, sep="\t", file=handle)


def main(argv: Optional[list[str]] = None) -> int:
    args = parse_args(argv)
    findings = check_reference(args.reference, args.grammar, args.majority)
    write_report(findings, args.report)
    message = summarise(findings, args.parameter, args.report)
    with open(
        args.summary, "w", encoding="utf-8", errors="surrogateescape"
    ) as handle:
        print(status(findings), file=handle)
        if message:
            print(message, file=handle)
    return 0


if __name__ == "__main__":
    sys.exit(main())
