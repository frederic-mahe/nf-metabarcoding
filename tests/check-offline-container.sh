#!/usr/bin/env bash
# [S83]: assert the air-gapped / pre-built-image container path resolves.
# A site whose compute nodes have no outbound network cannot use Seqera
# Wave (which builds the image at task start, [S08]). Instead it enables
# the engine and points `process.container` at a pre-pulled image in its
# `-c site.config` ([S75]), composing with the executor profile
# (`-profile slurm`) rather than an engine profile — so Wave is never
# turned on. This guards (a) that such a `-c` resolves the engine + the
# pre-built `process.container` with no `wave.enabled`, and (b) that a
# plain `-profile slurm` (no such `-c`) sets no `process.container`, so
# the default is unchanged. Config resolution only — offline execution
# stays a manual cluster smoke test, as for the rest of [S08]. Requires
# `nextflow` on PATH.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
readonly REPO_ROOT
cd "$REPO_ROOT"

readonly IMAGE='/scratch/shared/nf-metabarcoding/nf-metabarcoding.sif'

site_cfg="$(mktemp)"
readonly site_cfg
trap 'rm -f "${site_cfg}"' EXIT

# a site's air-gapped config: enable the engine + a pre-pulled image, no
# Wave. Composed with -profile slurm (executor only).
cat > "${site_cfg}" <<CFG
singularity.enabled    = true
singularity.autoMounts = true
process.container      = '${IMAGE}'
CFG

fail=0

# resolve the flattened config for a profile (optionally with the site
# `-c`) and assert the needle is present / absent.
#   assert <present|absent> <desc> <needle> <profile> [site-config]
assert() {
    local mode="$1"
    local desc="$2"
    local needle="$3"
    local profile="$4"
    local cfg="${5:-}"
    local -a cmd=(nextflow)
    if [ -n "${cfg}" ]; then
        cmd+=(-c "${cfg}")
    fi
    cmd+=(config main.nf -flat -profile "${profile}")
    local flat
    if ! flat="$("${cmd[@]}" 2>&1)"; then
        echo "FAIL: ${desc}: config did not resolve"
        echo "${flat}"
        fail=1
        return
    fi
    if [ "${mode}" = present ]; then
        if ! grep -qF -- "${needle}" <<<"${flat}"; then
            echo "FAIL: ${desc}: missing '${needle}'"
            fail=1
            return
        fi
    else
        if grep -qF -- "${needle}" <<<"${flat}"; then
            echo "FAIL: ${desc}: unexpected '${needle}'"
            fail=1
            return
        fi
    fi
    echo "OK: ${desc}"
}

# the air-gapped -c resolves the engine + the pre-built image, no Wave.
assert present "site -c sets the pre-built process.container" \
    "process.container = '${IMAGE}'" "slurm" "${site_cfg}"
assert present "site -c enables the singularity engine" \
    "singularity.enabled = true" "slurm" "${site_cfg}"
assert absent  "Wave stays off on the air-gapped path" \
    "wave.enabled = true" "slurm" "${site_cfg}"

# default -profile slurm (no site -c) sets no container — behaviour
# unchanged for everyone not on the air-gapped path.
assert absent  "plain -profile slurm sets no process.container" \
    "process.container" "slurm"

if [ "${fail}" -ne 0 ]; then
    echo "offline-container: FAILED"
    exit 1
fi
echo "offline-container: OK"
