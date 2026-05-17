#!/usr/bin/env bats
#
# Smoke tests for bin/discover_fastq.py invoked as a CLI. Unit-level
# coverage of the pattern table lives in tests/python/test_discover_fastq.py;
# this file only checks that the CLI wiring (argument parsing, TSV
# output, multi-folder handling, --extra-pattern flag) works end-to-end.
#
# COVERAGE: [S10], [S11], [S12]
#
# Run from the repository root:
#   bats tests/bin/

setup() {
    REPO_ROOT="$(git rev-parse --show-toplevel)"
    SCRIPT="${REPO_ROOT}/bin/discover_fastq.py"
    WORK="$(mktemp -d)"
}

teardown() {
    rm -rf "${WORK}"
}

@test "emits TSV (sample_id, r1, r2) for a canonical paired-end pair" {
    : > "${WORK}/Foo_L001_R1_001.fastq.gz"
    : > "${WORK}/Foo_L001_R2_001.fastq.gz"

    run python3 "${SCRIPT}" "${WORK}"

    [ "${status}" -eq 0 ]
    [ "$(echo "${output}" | wc -l)" -eq 1 ]
    expected_r1="${WORK}/Foo_L001_R1_001.fastq.gz"
    expected_r2="${WORK}/Foo_L001_R2_001.fastq.gz"
    [ "${output}" = "Foo	${expected_r1}	${expected_r2}" ]
}

@test "leaves the R2 column empty for a single-end fastq" {
    : > "${WORK}/loose_sample.fastq.gz"

    run python3 "${SCRIPT}" "${WORK}"

    [ "${status}" -eq 0 ]
    # trailing tab, then empty R2 column
    [ "${output}" = "loose_sample	${WORK}/loose_sample.fastq.gz	" ]
}

@test "walks multiple folders given on the command line" {
    mkdir -p "${WORK}/a" "${WORK}/b"
    : > "${WORK}/a/sampleA_R1.fastq.gz"
    : > "${WORK}/a/sampleA_R2.fastq.gz"
    : > "${WORK}/b/loose.fq"

    run python3 "${SCRIPT}" "${WORK}/a" "${WORK}/b"

    [ "${status}" -eq 0 ]
    [ "$(echo "${output}" | wc -l)" -eq 2 ]
    echo "${output}" | grep -q "^sampleA	"
    echo "${output}" | grep -q "^loose	"
}

@test "--extra-pattern overrides the canonical table" {
    : > "${WORK}/study42-mateA.fastq.gz"
    : > "${WORK}/study42-mateB.fastq.gz"

    run python3 "${SCRIPT}" --extra-pattern '*-mate{A,B}.fastq.gz' "${WORK}"

    [ "${status}" -eq 0 ]
    [ "$(echo "${output}" | wc -l)" -eq 1 ]
    echo "${output}" | grep -q "^study42	"
}

@test "warns on stderr when an R1 has no R2 partner" {
    : > "${WORK}/orphan_L001_R1_001.fastq.gz"

    run python3 "${SCRIPT}" "${WORK}"

    [ "${status}" -eq 0 ]
    # stdout: TSV with empty R2
    echo "${output}" | grep -q "orphan_L001_R1_001	"
    # stderr is captured into $output too because `run` merges 2>&1 by
    # default — assert the warning text appears somewhere.
    echo "${output}" | grep -q "treating as single-end"
}
