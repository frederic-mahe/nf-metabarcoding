#!/usr/bin/env bash
# [S76]: assert the self-contained `demo` profile is wired to the
# committed demo dataset so `nextflow run main.nf -profile demo` (or
# `-profile demo,singularity`) runs A -> B -> C out of the box, with no
# fixture generation and no required flags. Config resolution + asset
# presence only; the end-to-end run itself is exercised by the Part
# A->B->C nf-tests in tests/main.nf.test (same primers / shape) and by a
# manual `-profile demo` smoke run. Requires `nextflow` on PATH.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
readonly REPO_ROOT
cd "${REPO_ROOT}"

fail=0

assert_file() {
    local -r f="$1"
    if [ ! -s "${f}" ]; then
        echo "FAIL: bundled demo asset missing or empty: ${f}"
        fail=1
        return
    fi
    echo "OK: ${f} present"
}

assert_resolves() {
    local -r needle="$1"
    local flat
    if ! flat="$(nextflow config main.nf -flat -profile demo 2>&1)"; then
        echo "FAIL: 'nextflow config -profile demo' did not resolve"
        echo "${flat}"
        fail=1
        return
    fi
    if ! grep -qE -- "${needle}" <<<"${flat}"; then
        echo "FAIL: -profile demo is missing /${needle}/"
        fail=1
        return
    fi
    echo "OK: -profile demo -> /${needle}/"
}

# the committed dataset that makes the profile self-contained.
assert_file "assets/demo/demo_1.fastq"
assert_file "assets/demo/demo_2.fastq"
assert_file "assets/demo/reference.fasta"

# the profile points every required parameter at that dataset, so no
# flags are needed to run all three parts.
assert_resolves "params.fastq_folder = .*assets/demo"
assert_resolves "params.forward_primer = 'CCAGCASCYGCGGTAATTCC'"
assert_resolves "params.reverse_primer = 'ACTTTCGTTCTTGATYRA'"
assert_resolves "params.project_name = 'demo'"
assert_resolves "params.reference_dataset = .*assets/demo/reference.fasta"

# composes with a container engine for a one-command environment check.
if ! nextflow config main.nf -flat -profile demo,singularity >/dev/null 2>&1; then
    echo "FAIL: -profile demo,singularity did not resolve"
    fail=1
else
    echo "OK: -profile demo,singularity resolves"
fi

if [ "${fail}" -ne 0 ]; then
    echo "demo-profile: FAILED"
    exit 1
fi
echo "demo-profile: OK"
