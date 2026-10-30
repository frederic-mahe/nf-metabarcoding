#!/usr/bin/env bash
# [S110]: assert main.nf compiles on the Nextflow found on PATH. The
# --help entry point parses and compiles the whole script but schedules
# no process, so no tool, fixture or container is needed. Run by the CI
# `config-checks` matrix on a pinned 25.10.x and on the latest stable,
# whose strict parser (26.10+) rejects `return` in a workflow body.
# Requires `nextflow` on PATH.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
readonly REPO_ROOT
cd "${REPO_ROOT}"

fail=0

out=""
if ! out="$(nextflow run main.nf --help 2>&1)"; then
    echo "FAIL: 'nextflow run main.nf --help' exited non-zero"
    fail=1
fi
if grep -q "Script compilation failed" <<<"${out}"; then
    echo "FAIL: main.nf does not compile"
    fail=1
fi
if ! grep -q "nf-metabarcoding — swarm-based metabarcoding pipeline" <<<"${out}"; then
    echo "FAIL: --help did not print the usage banner"
    fail=1
fi

if [ "${fail}" -ne 0 ]; then
    echo "${out}"
    echo "compile: FAILED"
    exit 1
fi
echo "compile: OK"
