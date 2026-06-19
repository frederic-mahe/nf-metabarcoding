"""Release-version consistency check ([S77]).

The release version is declared in two files that must stay in sync:
``manifest.version`` in ``nextflow.config`` and ``version`` in
``CITATION.cff``. A release that bumps one but not the other would ship
citation metadata disagreeing with the workflow manifest. This check
asserts the two are identical (the automated guard behind the README
"Releasing" instructions).
"""

# COVERAGE: [S77]

from __future__ import annotations

import re
from pathlib import Path


def _manifest_version(repo_root: Path) -> str:
    text = (repo_root / "nextflow.config").read_text()
    # the manifest block's `version = '...'` — not the sibling
    # `nextflowVersion`, and not any `version` inside other scopes.
    block = re.search(r"manifest\s*\{(.*?)\}", text, re.S)
    assert block, "no manifest { } block in nextflow.config"
    match = re.search(
        r"^\s*version\s*=\s*['\"]([^'\"]+)['\"]", block.group(1), re.M
    )
    assert match, "no manifest.version in nextflow.config"
    return match.group(1).strip()


def _citation_version(repo_root: Path) -> str:
    text = (repo_root / "CITATION.cff").read_text()
    match = re.search(r"^version:\s*['\"]?([^'\"\s]+)['\"]?\s*$", text, re.M)
    assert match, "no version: in CITATION.cff"
    return match.group(1).strip()


def test_manifest_and_citation_versions_agree(repo_root: Path) -> None:
    # COVERAGE: [S77]
    manifest = _manifest_version(repo_root)
    citation = _citation_version(repo_root)
    assert manifest == citation, (
        f"nextflow.config manifest.version ({manifest!r}) and CITATION.cff "
        f"version ({citation!r}) must match — bump both (see the README "
        f"'Releasing' section)"
    )
