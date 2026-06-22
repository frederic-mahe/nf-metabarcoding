#!/bin/bash
#
# Build the per-sample read-count summary table ([S86]) — a port of the
# legacy genotoul read-tracking summary.
#
# Usage:
#   build_read_counts.sh <sample_ids_file> <output_tsv>
#
# Each sample's published Part A logs ([S19]) are read from the current
# working directory:
#   <id>_merging.log            (vsearch --fastq_mergepairs)
#   <id>_trimming_forward.log   (cutadapt forward pass)
#   <id>_trimming_reverse.log   (cutadapt reverse pass)
#
# The output is a tab-separated table:
#   samples  reads  assembled  F  R  passing
# one row per sample id (sorted), followed by a Total row. A count whose
# source log or source line is absent is recorded as 0 — every cell is
# numeric. Thousands separators are stripped.

set -euo pipefail

readonly IDS_FILE="${1:?usage: build_read_counts.sh <sample_ids_file> <output_tsv>}"
readonly OUTPUT="${2:?usage: build_read_counts.sh <sample_ids_file> <output_tsv>}"

# Extract a numeric count from a log file: the first line matching the
# awk pattern, taking the requested field, with thousands separators
# stripped. Prints 0 when the file is absent or no line matches.
extract() {
    local file=$1
    local pattern=$2
    local field=$3
    if [[ ! -f "${file}" ]]; then
        printf '0'
        return
    fi
    local value
    value="$(awk -v pat="${pattern}" -v fld="${field}" \
        '$0 ~ pat { gsub(/,/, "", $fld); print $fld; exit }' "${file}")"
    if [[ -z "${value}" ]]; then
        printf '0'
    else
        printf '%s' "${value}"
    fi
}

{
    printf 'samples\treads\tassembled\tF\tR\tpassing\n'

    total_reads=0
    total_assembled=0
    total_f=0
    total_r=0
    total_passing=0

    while IFS= read -r id; do
        [[ -z "${id}" ]] && continue

        local_reads="$(extract "${id}_merging.log" 'Pairs' 1)"
        local_assembled="$(extract "${id}_merging.log" 'Merged' 1)"
        local_f="$(extract "${id}_trimming_forward.log" 'Reads with adapters' 4)"
        local_r="$(extract "${id}_trimming_reverse.log" 'Reads with adapters' 4)"
        local_passing="$(extract "${id}_trimming_reverse.log" 'Reads written .passing filters.' 5)"

        printf '%s\t%d\t%d\t%d\t%d\t%d\n' \
            "${id}" \
            "${local_reads}" \
            "${local_assembled}" \
            "${local_f}" \
            "${local_r}" \
            "${local_passing}"

        total_reads=$(( total_reads + local_reads ))
        total_assembled=$(( total_assembled + local_assembled ))
        total_f=$(( total_f + local_f ))
        total_r=$(( total_r + local_r ))
        total_passing=$(( total_passing + local_passing ))
    done < <(sort "${IDS_FILE}")

    printf 'Total\t%d\t%d\t%d\t%d\t%d\n' \
        "${total_reads}" \
        "${total_assembled}" \
        "${total_f}" \
        "${total_r}" \
        "${total_passing}"
} > "${OUTPUT}"
