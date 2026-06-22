#!/usr/bin/env bash
# [S85]: assert the whole pipeline runs under `nextflow -stub-run` with
# NO bioinformatics tool present — i.e. every tool-invoking process
# declares a `stub:` that produces its declared outputs. This validates
# the Part A -> B -> C channel topology (joins, branches, scatter/gather)
# in seconds, without vsearch / swarm / cutadapt / mumu and without real
# data — the fast topology-CI / onboarding smoke check.
#
# Two layers:
#   1. static gate — every process module declares a `stub:`, except the
#      input-discovery / samplesheet-validation processes, which are pure
#      stdlib-Python glue that must run for real to bootstrap the
#      per-sample channel from the filesystem.
#   2. dynamic — run `-profile demo -stub-run` with the four tools
#      shadowed by stubs that `exit 1`, so if any process fell through to
#      its real script the run would fail. Python is left intact (the
#      discovery processes need it).
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
readonly REPO_ROOT
cd "$REPO_ROOT"

# processes that legitimately have no stub (real, tool-free discovery).
readonly -a STUB_EXEMPT=(
    discover_inputs
    discover_part_b_fasta
    validate_samplesheet
)

fail=0

# ---- layer 1: static stub-presence gate -------------------------------
is_exempt() {
    local name="$1"
    local e
    for e in "${STUB_EXEMPT[@]}"; do
        [ "${e}" = "${name}" ] && return 0
    done
    return 1
}

while IFS= read -r module; do
    name="$(sed -n 's/^process \([A-Za-z0-9_]*\).*/\1/p' "${module}" | head -n 1)"
    [ -z "${name}" ] && continue
    if is_exempt "${name}"; then
        continue
    fi
    if ! grep -qE '^[[:space:]]*stub:' "${module}"; then
        echo "FAIL: process '${name}' (${module#"${REPO_ROOT}"/}) declares no stub:"
        fail=1
    fi
done < <(grep -rl '^process ' "${REPO_ROOT}/modules/local" --include='*.nf')

if [ "${fail}" -eq 0 ]; then
    echo "OK: every non-discovery process declares a stub"
fi

# ---- layer 2: tool-free -stub-run of the demo topology ----------------
sandbox="$(mktemp -d)"
readonly sandbox
trap 'rm -rf "${sandbox}" "${REPO_ROOT}/demo_results" "${REPO_ROOT}/work" "${REPO_ROOT}"/.nextflow*' EXIT

# shadow the bioinformatics tools with stubs that fail loudly: if a real
# script runs, it calls one of these and the task dies.
for tool in vsearch swarm cutadapt mumu; do
    printf '#!/bin/sh\necho "FAIL: real %s was invoked under -stub-run" >&2\nexit 1\n' \
        "${tool}" > "${sandbox}/${tool}"
    chmod +x "${sandbox}/${tool}"
done

log="${sandbox}/stub-run.log"
if PATH="${sandbox}:${PATH}" nextflow run main.nf -profile demo -stub-run \
        > "${log}" 2>&1; then
    echo "OK: -profile demo -stub-run completed with no tool installed"
else
    echo "FAIL: -profile demo -stub-run did not complete"
    grep -iE 'FAIL: real|terminated|ERROR' "${log}" | head -n 5
    fail=1
fi

# the placeholder Part A/B/C artefacts must have been published.
for artefact in \
    demo_results/per_sample/demo.fas \
    demo_results/occurrence_table/demo_1_samples_table.tsv \
    demo_results/occurrence_table/demo_1_samples_table_assigned.tsv \
    demo_results/pipeline_info/software_versions.yml ; do
    if [ -e "${REPO_ROOT}/${artefact}" ]; then
        echo "OK: published ${artefact}"
    else
        echo "FAIL: stub-run did not publish ${artefact}"
        fail=1
    fi
done

if [ "${fail}" -ne 0 ]; then
    echo "stub-run: FAILED"
    exit 1
fi
echo "stub-run: OK"
