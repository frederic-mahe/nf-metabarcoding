#!/usr/bin/env bash
# [S07]: assert the `slurm` executor profile resolves. The slurm profile
# carries the only non-trivial config in the pipeline — per-process
# resource tiers with `task.attempt` escalation and dataset/reference
# size-scaled memory closures (some interpolating `as int` GB strings).
# A typo in one of those closures only surfaces on a real cluster
# otherwise, where the feedback loop is slow and expensive. This guards
# that (a) `nextflow config -profile slurm` parses every closure and
# resolves, (b) the slurm executor and the resourceLimits ceiling are
# wired, and (c) the profile composes with a dependency profile
# (`slurm,conda`). Requires `nextflow` on PATH; config resolution only —
# no scheduler and no job submission are exercised here (that stays a
# manual cluster smoke test, as for the rest of [S07] / [S08]).
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
readonly REPO_ROOT
cd "$REPO_ROOT"

fail=0

# resolve the flattened config for a profile and assert it contains the
# needle. A non-zero `nextflow config` exit means a closure failed to
# parse — reported as a failure with the resolver output.
#   assert_contains <desc> <needle> <profile>
assert_contains() {
    local desc="$1"
    local needle="$2"
    local profile="$3"
    local flat
    if ! flat="$(nextflow config main.nf -flat -profile "${profile}" 2>&1)"; then
        echo "FAIL: ${desc}: config did not resolve"
        echo "${flat}"
        fail=1
        return
    fi
    if ! grep -qF -- "${needle}" <<<"${flat}"; then
        echo "FAIL: ${desc}: missing '${needle}'"
        fail=1
        return
    fi
    echo "OK: ${desc}"
}

# the slurm profile wires the slurm executor and the resourceLimits cap.
assert_contains "slurm profile sets the slurm executor" \
    "process.executor = 'slurm'" "slurm"
assert_contains "slurm profile sets the resourceLimits ceiling" \
    "process.resourceLimits.cpus" "slurm"
assert_contains "slurm profile sets the submit rate limit" \
    "executor.submitRateLimit" "slurm"

# [S07]: slurm job arrays are an opt-in knob (--slurm_array_size), off by
# default so a plain `-profile slurm` keeps one job per task. null here
# means the `array` directive is omitted at runtime.
assert_contains "slurm profile leaves job arrays off by default" \
    "process.array = null" "slurm"

# the slurm profile composes with a dependency profile.
assert_contains "slurm composes with conda (executor)" \
    "process.executor = 'slurm'" "slurm,conda"
assert_contains "slurm composes with conda (conda enabled)" \
    "conda.enabled = true" "slurm,conda"

if [ "${fail}" -ne 0 ]; then
    echo "slurm-config: FAILED"
    exit 1
fi
echo "slurm-config: OK"
