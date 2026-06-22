#!/usr/bin/env bash
# [S87]: assert each vendored institutional cluster profile resolves to
# the expected wiring. Like check-slurm-config.sh / check-container-
# profiles.sh, this guards the profile *wiring* via `nextflow config` —
# a typo in a queue/clusterOptions closure or a wrong resourceLimits
# ceiling fails fast in CI instead of surfacing on a slow, expensive
# cluster run. It does NOT submit jobs: the partition/account *values*
# can only be confirmed on the hardware, so actual submission stays a
# manual cluster smoke test (per [S00], as for [S07] / [S08]).
#
# Each cluster profile implies slurm (it includes conf/slurm.config) and
# clamps resourceLimits to that cluster's largest node, then composes
# with a container engine (`-profile <cluster>,singularity`). Requires
# `nextflow` on PATH; config resolution only — no scheduler needed.
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

# Each cluster: profile name -> its resourceLimits memory ceiling (the
# largest node, as Nextflow renders it in `-flat` output).
declare -A CLUSTER_MEMORY=(
    [abims]='750 GB'
    [genotoul]='3.4 TB'
    [ifb_core]='252 GB'
    [meso]='4 TB'
    [saga]='6 TB'
)

for cluster in "${!CLUSTER_MEMORY[@]}"; do
    # implies slurm (includes conf/slurm.config) ...
    assert_contains "${cluster}: slurm executor" \
        "process.executor = 'slurm'" "${cluster}"
    # ... and clamps to that cluster's largest node ...
    assert_contains "${cluster}: resourceLimits memory ceiling" \
        "process.resourceLimits.memory = '${CLUSTER_MEMORY[$cluster]}'" "${cluster}"
    # ... and composes with a container engine.
    assert_contains "${cluster},singularity: engine enabled" \
        "singularity.enabled = true" "${cluster},singularity"
done

# meso derives partition + account per task from its own closures.
assert_contains "meso: per-task partition selector" \
    "task.memory > 1536.GB ? 'bigmem-cirad-dedicated' : 'cpu-dedicated'" "meso"
assert_contains "meso: per-task account routing" \
    "dedicated-cpu@cirad-normal" "meso"

# A plain run (no profile) must not silently enable the slurm executor.
flat="$(nextflow config main.nf -flat 2>&1)"
if grep -qF -- "process.executor = 'slurm'" <<<"${flat}"; then
    echo "FAIL: default (no profile) unexpectedly sets the slurm executor"
    fail=1
else
    echo "OK: default (no profile) does not set the slurm executor"
fi

if [ "${fail}" -ne 0 ]; then
    echo "cluster-profiles: FAILED"
    exit 1
fi
echo "cluster-profiles: OK"
