"""Repository-metadata consistency checks ([S80]).

The sibling of [S77]'s version-sync guard, but for the *links and
release notes* rather than the version string. Three pieces of metadata
drift independently and each ships a broken artefact when it does:

- ``manifest.homePage`` (``nextflow.config``) and ``url``
  (``CITATION.cff``) must name the same GitHub repository, or one of the
  two advertises a dead repo;
- ``nextflow_schema.json``'s ``$id`` URL must point at the repository's
  default branch (``manifest.defaultBranch``) and the same repo, or the
  published schema link 404s on a renamed/stale branch;
- ``CHANGELOG.md`` must exist and document the current
  ``manifest.version``, or a release ships with no changelog entry.

These are pure file-content assertions — no Nextflow run required.
"""

# COVERAGE: [S80]

from __future__ import annotations

import json
import re
from pathlib import Path


def _manifest_field(repo_root: Path, field: str) -> str:
    text = (repo_root / "nextflow.config").read_text()
    block = re.search(r"manifest\s*\{(.*?)\}", text, re.S)
    assert block, "no manifest { } block in nextflow.config"
    match = re.search(
        rf"^\s*{field}\s*=\s*['\"]([^'\"]+)['\"]", block.group(1), re.M
    )
    assert match, f"no manifest.{field} in nextflow.config"
    return match.group(1).strip()


def _citation_url(repo_root: Path) -> str:
    text = (repo_root / "CITATION.cff").read_text()
    match = re.search(r"^url:\s*['\"]?([^'\"\s]+)['\"]?\s*$", text, re.M)
    assert match, "no url: in CITATION.cff"
    return match.group(1).strip()


def _owner_repo(url: str) -> str:
    """Extract ``owner/repo`` from a github.com URL, ignoring scheme,
    a trailing slash, and a trailing ``.git``."""
    match = re.search(r"github\.com[:/]+([^/]+/[^/]+?)(?:\.git)?/?$", url)
    assert match, f"cannot parse owner/repo from URL: {url!r}"
    return match.group(1)


def _schema_id(repo_root: Path) -> str:
    schema = json.loads((repo_root / "nextflow_schema.json").read_text())
    sid = schema.get("$id")
    assert sid, "no $id in nextflow_schema.json"
    return sid


def test_homepage_and_citation_url_same_repo(repo_root: Path) -> None:
    # COVERAGE: [S80]
    home = _owner_repo(_manifest_field(repo_root, "homePage"))
    cite = _owner_repo(_citation_url(repo_root))
    assert home == cite, (
        f"manifest.homePage repo ({home!r}) and CITATION.cff url repo "
        f"({cite!r}) must name the same GitHub repository"
    )


def test_schema_id_points_at_default_branch(repo_root: Path) -> None:
    # COVERAGE: [S80]
    branch = _manifest_field(repo_root, "defaultBranch")
    assert branch, "manifest.defaultBranch must be set"
    repo = _owner_repo(_manifest_field(repo_root, "homePage"))
    sid = _schema_id(repo_root)
    needle = f"/{repo}/{branch}/"
    assert needle in sid, (
        f"nextflow_schema.json $id ({sid!r}) must reference the default "
        f"branch and repo (expected to contain {needle!r}) — a stale "
        f"branch name in $id ships a dead schema link"
    )


def test_changelog_documents_current_version(repo_root: Path) -> None:
    # COVERAGE: [S80]
    changelog = repo_root / "CHANGELOG.md"
    assert changelog.exists(), "CHANGELOG.md is missing"
    version = _manifest_field(repo_root, "version")
    assert version in changelog.read_text(), (
        f"CHANGELOG.md must carry an entry for the current "
        f"manifest.version ({version!r})"
    )
