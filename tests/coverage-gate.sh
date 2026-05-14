#!/usr/bin/env bash
#
# coverage-gate.sh
#
# Audits the SPECIFICATIONS.md <-> tests/COVERAGE.md <-> tests/**.nf.test
# triangle:
#   1. every [Sxx] ID declared in SPECIFICATIONS.md appears in
#      tests/COVERAGE.md
#   2. every [Sxx] ID mentioned in a tests/**.nf.test // COVERAGE:
#      comment is declared in SPECIFICATIONS.md
#
# Run from the repository root:
#   bash tests/coverage-gate.sh
#
# Exit status: 0 on success, 1 on any uncovered or unknown ID.

set -euo pipefail

REPO_ROOT="$(git rev-parse --show-toplevel)"
readonly REPO_ROOT
readonly SPEC_FILE="${REPO_ROOT}/SPECIFICATIONS.md"
readonly COVERAGE_FILE="${REPO_ROOT}/tests/COVERAGE.md"
readonly TESTS_DIR="${REPO_ROOT}/tests"

if [[ ! -r "${SPEC_FILE}" ]]; then
    echo "coverage-gate: cannot read ${SPEC_FILE}" >&2
    exit 1
fi
if [[ ! -r "${COVERAGE_FILE}" ]]; then
    echo "coverage-gate: cannot read ${COVERAGE_FILE}" >&2
    exit 1
fi

# IDs declared in SPECIFICATIONS.md (any [Sxx] in a backtick span)
declared_ids="$(grep --only-matching --extended-regexp '\[S[0-9]+\]' "${SPEC_FILE}" | sort --unique)"

# IDs referenced in COVERAGE.md
covered_ids="$(grep --only-matching --extended-regexp '\[S[0-9]+\]' "${COVERAGE_FILE}" | sort --unique)"

# IDs referenced from tests (nf-test, bats, and pytest sources)
test_ids="$(grep --recursive --no-filename --only-matching --extended-regexp \
    --include='*.nf.test' --include='*.bats' --include='test_*.py' \
    '\[S[0-9]+\]' "${TESTS_DIR}" 2>/dev/null | sort --unique || true)"

status=0

missing_in_coverage="$(comm -23 <(printf '%s\n' "${declared_ids}") <(printf '%s\n' "${covered_ids}"))"
if [[ -n "${missing_in_coverage}" ]]; then
    echo "coverage-gate: the following [Sxx] IDs are declared in SPECIFICATIONS.md"
    echo "               but missing from tests/COVERAGE.md:"
    awk '{ print "  " $0 }' <<< "${missing_in_coverage}"
    status=1
fi

unknown_in_tests="$(comm -23 <(printf '%s\n' "${test_ids}") <(printf '%s\n' "${declared_ids}"))"
if [[ -n "${unknown_in_tests}" ]]; then
    echo "coverage-gate: the following [Sxx] IDs are referenced from tests/"
    echo "               but not declared in SPECIFICATIONS.md:"
    awk '{ print "  " $0 }' <<< "${unknown_in_tests}"
    status=1
fi

unknown_in_coverage="$(comm -23 <(printf '%s\n' "${covered_ids}") <(printf '%s\n' "${declared_ids}"))"
if [[ -n "${unknown_in_coverage}" ]]; then
    echo "coverage-gate: the following [Sxx] IDs are referenced from tests/COVERAGE.md"
    echo "               but not declared in SPECIFICATIONS.md:"
    awk '{ print "  " $0 }' <<< "${unknown_in_coverage}"
    status=1
fi

if [[ "${status}" -eq 0 ]]; then
    declared_count="$(printf '%s\n' "${declared_ids}" | wc --lines | tr --delete ' ')"
    echo "coverage-gate: OK (${declared_count} [Sxx] IDs, all mapped)"
fi

exit "${status}"
