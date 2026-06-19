#!/usr/bin/env bash
# [S08]: assert each container engine profile resolves to the expected
# directives. Container *execution* (actually running a process inside a
# container) stays a manual cluster smoke test — it is not exercised
# here. This guards the profile *wiring* so a typo or a deleted directive
# fails fast in CI instead of surfacing on a cluster run. Requires
# `nextflow` on PATH; no bioinformatics tools needed (config resolution
# only).
#
# The four engine profiles (docker / podman / singularity / apptainer)
# each enable their own engine plus Wave, and point Wave at
# environment.yml (the single pinned source of truth, [S69]) so the
# image is built on the fly — no Dockerfile, no registry (DECISIONS.md
# D10).
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
readonly REPO_ROOT
cd "$REPO_ROOT"

fail=0

# assert that `nextflow config -flat -profile <profile>` resolves and
# contains the dotted `<key> = <value>` line.
assert_resolves() {
    local profile="$1"
    local needle="$2"
    local flat
    if ! flat="$(nextflow config main.nf -flat -profile "$profile" 2>&1)"; then
        echo "FAIL: 'nextflow config -profile ${profile}' did not resolve:"
        echo "${flat}"
        fail=1
        return
    fi
    if ! grep -qF -- "${needle}" <<<"${flat}"; then
        echo "FAIL: -profile ${profile} is missing '${needle}'"
        fail=1
        return
    fi
    echo "OK: -profile ${profile} -> ${needle}"
}

# assert a dotted key is absent (used to confirm the default profile
# does not silently enable containers / conda).
assert_absent() {
    local profile="$1"
    local needle="$2"
    local flat
    flat="$(nextflow config main.nf -flat ${profile:+-profile "$profile"} 2>&1)"
    if grep -qF -- "${needle}" <<<"${flat}"; then
        echo "FAIL: profile '${profile:-<none>}' unexpectedly sets '${needle}'"
        fail=1
        return
    fi
    echo "OK: profile '${profile:-<none>}' does not set '${needle}'"
}

# Each engine on. All four route the build through Wave + the pinned
# conda environment.
for engine in docker podman singularity apptainer; do
    assert_resolves "${engine}" "${engine}.enabled = true"
    assert_resolves "${engine}" "wave.enabled = true"
    assert_resolves "${engine}" "conda.enabled = true"
    assert_resolves "${engine}" "environment.yml"
done

# Singularity / apptainer additionally need autoMounts so host paths
# (work dir, inputs) are visible inside the container.
assert_resolves "singularity" "singularity.autoMounts = true"
assert_resolves "apptainer"   "apptainer.autoMounts = true"

# Composes with the slurm executor profile (the documented HPC entry
# point: `-profile slurm,singularity`).
assert_resolves "slurm,singularity" "process.executor = 'slurm'"
assert_resolves "slurm,singularity" "singularity.enabled = true"

# The default (no profile) must stay bare-PATH: no container engine, no
# Wave, no conda env creation.
assert_absent "" "wave.enabled = true"
assert_absent "" "conda.enabled = true"

if [ "${fail}" -ne 0 ]; then
    echo "container-profiles: FAILED"
    exit 1
fi
echo "container-profiles: OK"
