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
    [abims]='2 TB'
    [genotoul]='3.9 TB'
    [ifb_core]='2 TB'
    [meso]='3.7 TB'
    [saga]='5.9 TB'
)

# The dependency engine each cluster composes with. The needle is the
# engine's `enabled` flag in the resolved config. Genotoul offers
# apptainer/singularity as MODULES (loaded by its beforeScript, asserted
# separately below); the rest expose singularity on PATH.
declare -A CLUSTER_ENGINE=(
    [abims]='singularity'
    [genotoul]='apptainer'
    [ifb_core]='singularity'
    [meso]='singularity'
    [saga]='singularity'
)

for cluster in "${!CLUSTER_MEMORY[@]}"; do
    engine="${CLUSTER_ENGINE[$cluster]}"
    # implies slurm (includes conf/slurm.config) ...
    assert_contains "${cluster}: slurm executor" \
        "process.executor = 'slurm'" "${cluster}"
    # ... and clamps to that cluster's largest node ...
    assert_contains "${cluster}: resourceLimits memory ceiling" \
        "process.resourceLimits.memory = '${CLUSTER_MEMORY[$cluster]}'" "${cluster}"
    # ... composes with its dependency engine ...
    assert_contains "${cluster},${engine}: engine enabled" \
        "${engine}.enabled = true" "${cluster},${engine}"
    # ... and defaults the per-sample fan-out to slurm job arrays ([S87]).
    assert_contains "${cluster}: job-array size default" \
        "process.array = 50" "${cluster}"
done

# abims routes memory-heavy tasks to its bigmem partition (which has more
# memory than fast/long), otherwise picks fast/long by wall-time.
assert_contains "abims: bigmem memory routing" \
    "task.memory > 1400.GB ? 'bigmem'" "abims"

# genotoul loads its container engine from a module via beforeScript.
assert_contains "genotoul: apptainer module beforeScript" \
    "module load containers/Apptainer/1.4.1" "genotoul,apptainer"

# meso derives partition + account per task from its own closures.
assert_contains "meso: per-task partition selector" \
    "task.memory > 1400.GB ? 'bigmem-cirad-dedicated' : 'cpu-dedicated'" "meso"
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
