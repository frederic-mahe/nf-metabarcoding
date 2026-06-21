#!/usr/bin/env bash
# [S82]: assert `cleanup` defaults to `false`. Retaining per-task work/
# directories on success is what makes `-resume` work across separate
# invocations (the large-dataset workflow) and keeps a run inspectable;
# auto-clean is opt-in via a `-c site.config`. A regression to
# `cleanup = true` would silently break cross-invocation `-resume`, so
# guard the resolved default. Config resolution only — Nextflow's actual
# cleanup behaviour is upstream and not re-tested ([S00]). Requires
# `nextflow` on PATH.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
readonly REPO_ROOT
cd "$REPO_ROOT"

fail=0

# assert the resolved (flattened) config for a profile contains the
# needle.  assert_contains <desc> <needle> [profile]
assert_contains() {
    local desc="$1"
    local needle="$2"
    local profile="${3:-}"
    local -a cmd=(nextflow config main.nf -flat)
    if [ -n "${profile}" ]; then
        cmd+=(-profile "${profile}")
    fi
    local flat
    if ! flat="$("${cmd[@]}" 2>&1)"; then
        echo "FAIL: ${desc}: config did not resolve"
        echo "${flat}"
        fail=1
        return
    fi
    if ! grep -qF -- "${needle}" <<<"${flat}"; then
        echo "FAIL: ${desc}: missing '${needle}'"
        echo "${flat}" | grep -E '^cleanup' || echo "  (no cleanup line resolved)"
        fail=1
        return
    fi
    echo "OK: ${desc}"
}

# default (no profile) and the test profile must both resolve cleanup off.
assert_contains "cleanup defaults to false (no profile)" "cleanup = false"
assert_contains "cleanup stays false under -profile test" "cleanup = false" "test"

if [ "${fail}" -ne 0 ]; then
    echo "cleanup-default: FAILED"
    exit 1
fi
echo "cleanup-default: OK"
