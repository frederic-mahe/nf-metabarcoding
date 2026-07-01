#!/usr/bin/env bash
# [S07]: assert the slurm-only resource/account params are ALWAYS defined,
# including on a plain `nextflow run` with no slurm profile. main.nf reads
# params.slurm_account (the [S92] account check) and params.dataset_size_gb
# / params.reference_size_gb (the [S79] size warnings) unconditionally, on
# every entry point. Those params carry their runtime defaults in the
# `slurm` profile (conf/slurm.config), so unless they are also declared at
# top level a local run trips Nextflow's "Access to undefined parameter"
# warning for each — the same reason params.require_slurm_account is
# declared in nextflow.config's top-level params block. This guards that
# the top-level declarations stay in place: `nextflow config` with NO
# profile must resolve all three (to their `null` unset sentinel).
# Config resolution only — no scheduler, no job submission.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
readonly REPO_ROOT
cd "$REPO_ROOT"

fail=0

# resolve the flattened config with NO profile and assert the param is
# present (defined). A missing line means the param is only declared in
# the slurm profile and would warn on a local run.
#   assert_defined <param>
assert_defined() {
    local param="$1"
    local flat
    if ! flat="$(nextflow config main.nf -flat 2>/dev/null)"; then
        echo "FAIL: ${param}: config did not resolve"
        fail=1
        return
    fi
    if ! grep -qE -- "^params\.${param} = " <<<"${flat}"; then
        echo "FAIL: params.${param} is not defined without a slurm profile"
        fail=1
        return
    fi
    echo "OK: params.${param} is defined without a slurm profile"
}

assert_defined "slurm_account"
assert_defined "dataset_size_gb"
assert_defined "reference_size_gb"

if [ "${fail}" -ne 0 ]; then
    echo "local-params: FAILED"
    exit 1
fi
echo "local-params: OK"
