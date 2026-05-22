#!/usr/bin/env bats
#
# Unit tests for bin/reverse_complement.sh — the IUPAC-aware
# reverse-complement helper that trim_primers ([S01], Part A) calls
# to derive the reverse primer's revcomp before feeding it to
# cutadapt's --adapter search.
#
# COVERAGE: [S01]
#
# Run from the repository root:
#   bats tests/bin/

setup() {
    REPO_ROOT="$(git rev-parse --show-toplevel)"
    SCRIPT="${REPO_ROOT}/bin/reverse_complement.sh"
}

@test "reverse-complements a standard ACGT primer (argument form)" {
    result="$(bash "${SCRIPT}" CCAGCASCYGCGGTAATTCC)"
    [ "${result}" = "GGAATTACCGCRGSTGCTGG" ]
}

@test "reverse-complements via stdin when no argument is given" {
    result="$(echo "ACTTTCGTTCTTGATYRA" | bash "${SCRIPT}")"
    [ "${result}" = "TYRATCAAGAACGAAAGT" ]
}

@test "preserves case (lowercase in, lowercase out)" {
    result="$(bash "${SCRIPT}" acgt)"
    [ "${result}" = "acgt" ]
}

@test "passes N through (N is its own complement)" {
    result="$(bash "${SCRIPT}" ACGNT)"
    [ "${result}" = "ANCGT" ]
}

@test "handles full IUPAC alphabet" {
    # Complement: R<->Y, K<->M, S/W self-complement, B<->V, D<->H
    #   RYKMSWBVDH -> YRMKSWVBHD, then reversed -> DHBVWSKMRY
    result="$(bash "${SCRIPT}" RYKMSWBVDH)"
    [ "${result}" = "DHBVWSKMRY" ]
}

@test "empty input returns empty output" {
    result="$(printf '' | bash "${SCRIPT}")"
    [ -z "${result}" ]
}
