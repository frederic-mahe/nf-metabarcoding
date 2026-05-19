#!/usr/bin/env bats
#
# Smoke tests for bin/discover_fasta.py invoked as a CLI. Unit-level
# coverage of the discovery rules lives in
# tests/python/test_discover_fasta.py; this file only checks the CLI
# wiring (argument parsing, TSV output, multi-folder handling, the
# _notmerged skip).
#
# COVERAGE: [S27]
#
# Run from the repository root:
#   bats tests/bin/

setup() {
    REPO_ROOT="$(git rev-parse --show-toplevel)"
    SCRIPT="${REPO_ROOT}/bin/discover_fasta.py"
    WORK="$(mktemp -d)"
}

teardown() {
    rm -rf "${WORK}"
}

@test "emits TSV (sample_id, fasta) for every .fas in a folder" {
    printf '>seq\nACGT\n' > "${WORK}/A.fas"
    printf '>seq\nACGT\n' > "${WORK}/B.fas"

    run python3 "${SCRIPT}" "${WORK}"

    [ "${status}" -eq 0 ]
    [ "$(echo "${output}" | wc -l)" -eq 2 ]
    echo "${output}" | grep -q "^A	${WORK}/A.fas$"
    echo "${output}" | grep -q "^B	${WORK}/B.fas$"
}

@test "skips _notmerged.fas artefacts" {
    printf '>seq\nACGT\n' > "${WORK}/A.fas"
    printf '>seq\nACGT\n' > "${WORK}/A_notmerged.fas"

    run python3 "${SCRIPT}" "${WORK}"

    [ "${status}" -eq 0 ]
    [ "$(echo "${output}" | wc -l)" -eq 1 ]
    echo "${output}" | grep -q "^A	${WORK}/A.fas$"
}

@test "walks multiple folders given on the command line" {
    mkdir -p "${WORK}/a" "${WORK}/b"
    printf '>seq\nACGT\n' > "${WORK}/a/sampleA.fas"
    printf '>seq\nACGT\n' > "${WORK}/b/sampleB.fas"

    run python3 "${SCRIPT}" "${WORK}/a" "${WORK}/b"

    [ "${status}" -eq 0 ]
    [ "$(echo "${output}" | wc -l)" -eq 2 ]
    echo "${output}" | grep -q "^sampleA	"
    echo "${output}" | grep -q "^sampleB	"
}

@test "aborts on duplicate sample IDs and lists offending paths" {
    mkdir -p "${WORK}/run1" "${WORK}/run2"
    printf '>seq\nACGT\n' > "${WORK}/run1/A.fas"
    printf '>seq\nACGT\n' > "${WORK}/run2/A.fas"

    run python3 "${SCRIPT}" "${WORK}/run1" "${WORK}/run2"

    [ "${status}" -ne 0 ]
    echo "${output}" | grep -qi "duplicate"
    echo "${output}" | grep -q "${WORK}/run1/A.fas"
    echo "${output}" | grep -q "${WORK}/run2/A.fas"
}
