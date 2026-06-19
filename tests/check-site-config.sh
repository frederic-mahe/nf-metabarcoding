#!/usr/bin/env bash
# [S75]: assert the site-config override mechanism resolves. A site
# adapts the pipeline to its cluster by copying conf/site.config.example,
# editing it, and passing it with `-c my-site.config` (native Nextflow) —
# no pipeline edit required. This guards (a) that the shipped example
# template still parses and resolves, (b) that its overrides actually win
# over the slurm-profile defaults, and (c) that the slurm clusterOptions
# closure is wired to params.slurm_clusterOptions (the QoS / constraints
# passthrough). Requires `nextflow` on PATH; config resolution only — no
# scheduler and no job submission are exercised here (that stays a manual
# cluster smoke test, as for the rest of [S07] / [S08]).
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
readonly REPO_ROOT
cd "$REPO_ROOT"

readonly EXAMPLE="conf/site.config.example"

fail=0

# resolve the config for a given profile (and optional `-c` site file)
# and assert the flat output contains the needle. The `-c` flag is a
# global launcher option, so it must precede the `config` subcommand.
#   assert_contains <desc> <needle> <profile> [site-config-file]
assert_contains() {
    local desc="$1"
    local needle="$2"
    local profile="$3"
    local cfg="${4:-}"
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
    if ! grep -qF -- "${needle}" <<<"${flat}"; then
        echo "FAIL: ${desc}: missing '${needle}'"
        fail=1
        return
    fi
    echo "OK: ${desc}"
}

# the example template must exist and parse.
if [ ! -r "${EXAMPLE}" ]; then
    echo "FAIL: ${EXAMPLE} is missing"
    exit 1
fi

# default: the slurm profile exposes the passthrough param (null default)
# and the clusterOptions closure reads it.
assert_contains "slurm exposes slurm_clusterOptions (null default)" \
    "params.slurm_clusterOptions = null" "slurm"
assert_contains "clusterOptions closure reads slurm_clusterOptions" \
    "slurm_clusterOptions" "slurm"

# the example template resolves and its overrides win over the profile
# defaults when layered with `-c`.
assert_contains "site config overrides slurm_clusterOptions" \
    "params.slurm_clusterOptions = '--qos=long --constraint=haswell'" \
    "slurm" "${EXAMPLE}"
assert_contains "site config sets the singularity cache dir" \
    "singularity.cacheDir" \
    "slurm,singularity" "${EXAMPLE}"

if [ "${fail}" -ne 0 ]; then
    echo "site-config: FAILED"
    exit 1
fi
echo "site-config: OK"
