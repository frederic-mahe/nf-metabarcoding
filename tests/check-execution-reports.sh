#!/usr/bin/env bash
# [S84]: assert the per-run execution reports are wired. Nextflow's
# timeline / report / trace / dag must be enabled and target
# <outdir>/pipeline_info/ (the sibling of software_versions.yml), because
# the trace + report are how a site reads real per-step peak_rss /
# realtime to tune resources on a new cluster. This guards the config
# wiring — that the four reports are enabled and point at
# pipeline_info/. Nextflow actually generating the files (and following
# the run-time --outdir / -params-file) is upstream behaviour, verified
# by the demo / cluster smoke runs, not re-tested here ([S00]); note that
# `nextflow config` always renders the parse-time default output root
# ('results'), so this checks the default path. Requires `nextflow` on
# PATH.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
readonly REPO_ROOT
cd "$REPO_ROOT"

flat="$(nextflow config main.nf -flat 2>&1)" || {
    echo "FAIL: config did not resolve"
    echo "${flat}"
    exit 1
}

fail=0

assert_line() {
    local desc="$1"
    local needle="$2"
    if grep -qF -- "${needle}" <<<"${flat}"; then
        echo "OK: ${desc}"
    else
        echo "FAIL: ${desc}: missing '${needle}'"
        fail=1
    fi
}

# all four report types enabled
assert_line "timeline enabled" "timeline.enabled = true"
assert_line "report enabled"   "report.enabled = true"
assert_line "trace enabled"    "trace.enabled = true"
assert_line "dag enabled"      "dag.enabled = true"

# each report targets <outdir>/pipeline_info/ (default outdir = results)
assert_line "timeline → pipeline_info" "timeline.file = 'results/pipeline_info/execution_timeline.html'"
assert_line "report → pipeline_info"   "report.file = 'results/pipeline_info/execution_report.html'"
assert_line "trace → pipeline_info"    "trace.file = 'results/pipeline_info/execution_trace.txt'"
assert_line "dag → pipeline_info"      "dag.file = 'results/pipeline_info/pipeline_dag.html'"

# the trace carries the resource-tuning columns ([S84] use case)
assert_line "trace records peak_rss" "peak_rss"
assert_line "trace records peak_vmem" "peak_vmem"

if [ "${fail}" -ne 0 ]; then
    echo "execution-reports: FAILED"
    exit 1
fi
echo "execution-reports: OK"
