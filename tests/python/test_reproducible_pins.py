"""Repo-level reproducibility checks on dependency pins ([S69]).

A floating pin (``>=``) lets an upstream tool release silently change
output and break the byte-exact characterization tests
([S22]/[S35]/[S39]/[S44]). These checks assert that:

* ``environment.yml`` pins every conda dependency (including
  ``bioconda::mumu``) to an exact ``=`` version;
* ``environment.yml`` declares a ``bioconda::mumu`` package (mumu is
  now distributed on bioconda, so it is resolved by the ``conda``
  profile rather than built from source);
* the ``vsearch`` / ``swarm`` / ``cutadapt`` / ``mumu`` pins in
  ``environment.yml`` and the CI workflow are identical.
"""

# COVERAGE: [S69]

from __future__ import annotations

import re
from pathlib import Path

# tools whose versions must agree between environment.yml and CI.
SHARED_TOOLS = ("vsearch", "swarm", "cutadapt", "mumu")


def _dependencies_block(text: str) -> list[str]:
    """Return the ``- entry`` items under the ``dependencies:`` key."""
    entries: list[str] = []
    in_deps = False
    for line in text.splitlines():
        if re.match(r"^\s*dependencies:\s*$", line):
            in_deps = True
            continue
        if not in_deps:
            continue
        stripped = line.strip()
        # a non-indented, non-list line ends the block (next top key).
        if line and not line[0].isspace() and not stripped.startswith("-"):
            break
        if stripped.startswith("- "):
            entries.append(stripped[2:].strip())
    return entries


def _parse_dep(entry: str) -> tuple[str, str]:
    """Split ``channel::name=ver`` into ``(name, '=ver')``."""
    body = entry.split("::", 1)[-1]
    match = re.match(r"([A-Za-z0-9_.-]+)\s*(.*)$", body)
    assert match, f"unparseable dependency entry: {entry!r}"
    return match.group(1), match.group(2).strip()


def _env_deps(repo_root: Path) -> dict[str, str]:
    text = (repo_root / "environment.yml").read_text()
    return dict(_parse_dep(e) for e in _dependencies_block(text))


def _ci_pins(repo_root: Path) -> dict[str, str]:
    text = (repo_root / ".github/workflows/test.yml").read_text()
    pins: dict[str, str] = {}
    for tool in SHARED_TOOLS:
        match = re.search(rf"\b{tool}\s*([=<>!~]=?\s*[\w.]+)", text)
        if match:
            pins[tool] = match.group(1).replace(" ", "")
    return pins


def test_environment_pins_are_exact(repo_root: Path) -> None:
    deps = _env_deps(repo_root)
    assert deps, "no dependencies parsed from environment.yml"
    for name, spec in deps.items():
        assert re.fullmatch(r"=[\w.]+", spec), (
            f"{name} must use an exact '=' pin, got {spec!r}"
        )


def test_environment_pins_bioconda_mumu(repo_root: Path) -> None:
    entries = _dependencies_block(
        (repo_root / "environment.yml").read_text()
    )
    mumu = [e for e in entries if "mumu" in e]
    assert mumu, (
        "mumu is now distributed on bioconda and must be declared as a "
        "pinned conda dependency ([S69]) so the conda profile resolves it"
    )
    assert all(re.fullmatch(r"bioconda::mumu=[\w.]+", e) for e in mumu), (
        f"mumu must be pinned as 'bioconda::mumu=X.Y.Z', got {mumu}"
    )


def test_ci_and_env_pins_agree(repo_root: Path) -> None:
    env = _env_deps(repo_root)
    ci = _ci_pins(repo_root)
    for tool in SHARED_TOOLS:
        assert tool in env, f"{tool} missing from environment.yml"
        assert tool in ci, f"{tool} missing from CI create-args"
        assert ci[tool] == env[tool], (
            f"{tool}: CI pin {ci[tool]!r} != env pin {env[tool]!r}"
        )
