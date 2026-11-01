"""Every parameter the workflow reads is declared in nextflow.config."""

# COVERAGE: [S118]

import re
from pathlib import Path


def _strip_comments(text: str) -> str:
    """Drop `/* ... */` blocks and `//` line comments (not `://`)."""
    text = re.sub(r"/\*.*?\*/", "", text, flags=re.S)
    return re.sub(r"(?<!:)//.*", "", text)


def _declared(config: Path) -> set[str]:
    """Names assigned in the top-level `params { }` block."""
    lines = config.read_text().splitlines()
    start = lines.index("params {")
    end = lines.index("}", start)
    block = _strip_comments("\n".join(lines[start + 1:end]))
    return set(re.findall(r"^\s*(\w+)\s*=", block, flags=re.M))


def _read(repo_root: Path) -> dict[str, set[str]]:
    """`params.<name>` reads in the workflow sources, by name."""
    sources = [repo_root / "main.nf"]
    for folder in ("modules", "subworkflows"):
        sources += sorted((repo_root / folder).rglob("*.nf"))
    reads: dict[str, set[str]] = {}
    for source in sources:
        text = _strip_comments(source.read_text())
        for name in re.findall(r"\bparams\.(\w+)", text):
            reads.setdefault(name, set()).add(
                str(source.relative_to(repo_root))
            )
    return reads


def test_every_param_read_is_declared(repo_root: Path) -> None:
    # [S118]: an undeclared param read logs `Access to undefined
    # parameter`; declare it (null when it has no default) instead.
    declared = _declared(repo_root / "nextflow.config")
    missing = {
        name: sorted(files)
        for name, files in _read(repo_root).items()
        if name not in declared
    }
    assert missing == {}


def test_input_mode_selectors_are_declared(repo_root: Path) -> None:
    # [S118]: validate_params() reads every selector on every run.
    declared = _declared(repo_root / "nextflow.config")
    assert {
        "occurrence_table", "representatives_fasta", "input",
        "fasta_folder", "fastq_folder", "accession",
    } <= declared
