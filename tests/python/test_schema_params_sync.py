"""Drift guard: nextflow_schema.json <-> nextflow.config params.

Phase 1 of the nf-schema migration adds a ``nextflow_schema.json`` that
mirrors the parameter surface declared in ``nextflow.config`` — both the
defaults in the top-level ``params { }`` block and the params injected by
the profiles (slurm / modules / demo). nf-schema does NOT push the
schema's defaults back into ``params`` (defaults still live in the
config), so the two files can silently drift. This check keeps them
honest, in the spirit of the manifest <-> CITATION version check:

  1. every key declared in the ``params { }`` block has a matching
     schema property whose ``default`` equals the config default
     (a ``= null`` config default maps to "no ``default`` key");
  2. every schema property is a real parameter the workflow reads
     (``params.<name>`` appears somewhere in the ``*.nf`` / ``*.config``
     sources), so the schema never accumulates dead entries.

Parameters that are referenced/profile-injected but have no top-level
default (``fastq_folder``, ``forward_primer``, ``reverse_primer``, the
slurm/modules knobs) are covered by check (2) only — they legitimately
have no entry in the base ``params { }`` block.
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any


def _params_block(config_text: str) -> str:
    """Return the body of the top-level ``params { }`` block.

    The closing brace of the block sits in column 0; values inside (e.g.
    ``fastq_pattern = "*_1_{1,2}.fastq.gz"``) contain their own braces,
    so a non-greedy ``{...}`` match would stop too early. Slice from the
    ``params {`` line to the first line that is a bare ``}``.
    """
    lines = config_text.splitlines()
    start = next(
        (i for i, ln in enumerate(lines) if re.match(r"^params\s*\{", ln)),
        None,
    )
    assert start is not None, "no `params {` block in nextflow.config"
    end = next(
        (i for i in range(start + 1, len(lines))
         if re.match(r"^\}", lines[i])),
        None,
    )
    assert end is not None, "unterminated `params {` block in nextflow.config"
    return "\n".join(lines[start + 1:end])


def _parse_value(raw: str) -> Any:
    """Parse a Groovy scalar literal into its Python equivalent."""
    raw = raw.strip()
    if raw and raw[0] in "'\"":
        quote = raw[0]
        return raw[1:raw.index(quote, 1)]
    # strip a trailing inline comment on unquoted values
    raw = raw.split("//", 1)[0].strip()
    if raw == "true":
        return True
    if raw == "false":
        return False
    if raw == "null":
        return None
    try:
        return int(raw)
    except ValueError:
        return float(raw)


def _config_param_defaults(config_text: str) -> dict[str, Any]:
    body = _params_block(config_text)
    defaults: dict[str, Any] = {}
    for line in body.splitlines():
        match = re.match(r"^\s*([A-Za-z_]\w*)\s*=\s*(.+?)\s*$", line)
        if not match:
            continue  # comment / blank line
        defaults[match.group(1)] = _parse_value(match.group(2))
    return defaults


def _schema_properties(schema: dict[str, Any]) -> dict[str, dict[str, Any]]:
    props: dict[str, dict[str, Any]] = {}
    props.update(schema.get("properties", {}))
    for group in schema.get("$defs", {}).values():
        props.update(group.get("properties", {}))
    return props


def _referenced_params(repo_root: Path) -> set[str]:
    names: set[str] = set()
    pattern = re.compile(r"params\.([A-Za-z_]\w*)")
    for path in [*repo_root.rglob("*.nf"), *repo_root.rglob("*.config")]:
        if ".nf-test" in path.parts:
            continue
        names.update(pattern.findall(path.read_text()))
    return names


def test_config_defaults_match_schema(repo_root: Path) -> None:
    schema = json.loads((repo_root / "nextflow_schema.json").read_text())
    props = _schema_properties(schema)
    config_text = (repo_root / "nextflow.config").read_text()
    defaults = _config_param_defaults(config_text)

    for name, value in defaults.items():
        assert name in props, (
            f"params.{name} is declared in nextflow.config but missing from "
            f"nextflow_schema.json"
        )
        prop = props[name]
        if value is None:
            assert "default" not in prop, (
                f"params.{name} defaults to null in nextflow.config but the "
                f"schema declares default {prop.get('default')!r}"
            )
        else:
            assert prop.get("default") == value, (
                f"params.{name} default mismatch: nextflow.config has "
                f"{value!r}, nextflow_schema.json has {prop.get('default')!r}"
            )


def test_schema_properties_are_used(repo_root: Path) -> None:
    schema = json.loads((repo_root / "nextflow_schema.json").read_text())
    props = _schema_properties(schema)
    referenced = _referenced_params(repo_root)

    for name in props:
        assert name in referenced, (
            f"nextflow_schema.json declares '{name}' but no `params.{name}` "
            f"is read anywhere in the *.nf / *.config sources (dead entry?)"
        )
