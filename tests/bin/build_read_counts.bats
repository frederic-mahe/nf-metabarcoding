#!/usr/bin/env bats
#
# Unit tests for bin/build_read_counts.sh — the per-sample read-count
# summary table builder ([S86], a port of the legacy genotoul
# read-tracking summary).
#
# The script reads each sample's published Part A logs from the current
# working directory (<id>_merging.log, <id>_trimming_forward.log,
# <id>_trimming_reverse.log) and emits a TSV:
#
#   samples  reads  assembled  F  R  passing
#
# one row per sample id listed in the ids file, plus a trailing Total
# row. Absent logs / missing lines count as 0.
#
# COVERAGE: [S86]
#
# Run from the repository root:
#   bats tests/bin/

setup() {
    REPO_ROOT="$(git rev-parse --show-toplevel)"
    SCRIPT="${REPO_ROOT}/bin/build_read_counts.sh"
    WORK="$(mktemp -d)"
    cd "${WORK}"
}

teardown() {
    rm -rf "${WORK}"
}

# Write a vsearch-style merging log with the given Pairs / Merged counts.
write_merging() {
    local id=$1
    local pairs=$2
    local merged=$3
    {
        printf 'vsearch v2.31.0\n'
        printf 'vsearch --fastq_mergepairs %s_1.fastq.gz --reverse ...\n' "${id}"
        printf 'Started  2026-06-18\n'
        printf '   %s  Pairs\n' "${pairs}"
        printf '   %s  Merged (100.0%%)\n' "${merged}"
        printf '   0  Not merged (0.0%%)\n'
    } > "${id}_merging.log"
}

# Write a cutadapt-style trimming log. $2 = "Reads with adapters" count,
# $3 = "Reads written (passing filters)" count.
write_trimming() {
    local file=$1
    local adapters=$2
    local passing=$3
    {
        printf 'This is cutadapt 5.2\n'
        printf '=== Summary ===\n'
        printf 'Total reads processed:                       %s\n' "${adapters}"
        printf 'Reads with adapters:                         %s (100.0%%)\n' "${adapters}"
        printf '== Read fate breakdown ==\n'
        printf 'Reads written (passing filters):             %s (100.0%%)\n' "${passing}"
    } > "${file}"
}

@test "emits the canonical header" {
    printf 'A\n' > ids.txt
    write_merging A 10 9
    write_trimming A_trimming_forward.log 9 9
    write_trimming A_trimming_reverse.log 5 4
    run bash "${SCRIPT}" ids.txt out.tsv
    [ "${status}" -eq 0 ]
    run head -n 1 out.tsv
    [ "${output}" = "$(printf 'samples\treads\tassembled\tF\tR\tpassing')" ]
}

@test "extracts the five counts for one fully-populated sample" {
    printf 'DIV_T_PCR1a_S253\n' > ids.txt
    write_merging DIV_T_PCR1a_S253 3849 3758
    write_trimming DIV_T_PCR1a_S253_trimming_forward.log 3748 3748
    write_trimming DIV_T_PCR1a_S253_trimming_reverse.log 2172 1935
    run bash "${SCRIPT}" ids.txt out.tsv
    [ "${status}" -eq 0 ]
    run grep -P '^DIV_T_PCR1a_S253\t' out.tsv
    [ "${output}" = "$(printf 'DIV_T_PCR1a_S253\t3849\t3758\t3748\t2172\t1935')" ]
}

@test "absent merging log and missing trimming summary count as zero" {
    # An empty sample: merging log has 0 Pairs / 0 Merged and cutadapt
    # produced no Summary block (no adapter / passing lines).
    printf 'EMPTY\n' > ids.txt
    write_merging EMPTY 0 0
    printf 'This is cutadapt 5.2\nNo reads processed!\n' > EMPTY_trimming_forward.log
    printf 'This is cutadapt 5.2\nNo reads processed!\n' > EMPTY_trimming_reverse.log
    run bash "${SCRIPT}" ids.txt out.tsv
    [ "${status}" -eq 0 ]
    run grep -P '^EMPTY\t' out.tsv
    [ "${output}" = "$(printf 'EMPTY\t0\t0\t0\t0\t0')" ]
}

@test "a single-end sample with no merging log gets zero reads/assembled" {
    printf 'SE\n' > ids.txt
    # no SE_merging.log
    write_trimming SE_trimming_forward.log 100 100
    write_trimming SE_trimming_reverse.log 80 70
    run bash "${SCRIPT}" ids.txt out.tsv
    [ "${status}" -eq 0 ]
    run grep -P '^SE\t' out.tsv
    [ "${output}" = "$(printf 'SE\t0\t0\t100\t80\t70')" ]
}

@test "strips thousands separators from tool output" {
    printf 'BIG\n' > ids.txt
    write_merging BIG 124,333 122,415
    write_trimming BIG_trimming_forward.log 122,143 122,143
    write_trimming BIG_trimming_reverse.log 35,328 26,669
    run bash "${SCRIPT}" ids.txt out.tsv
    [ "${status}" -eq 0 ]
    run grep -P '^BIG\t' out.tsv
    [ "${output}" = "$(printf 'BIG\t124333\t122415\t122143\t35328\t26669')" ]
}

@test "appends a Total row summing every numeric column" {
    printf 'A\nB\n' > ids.txt
    write_merging A 10 9
    write_trimming A_trimming_forward.log 9 9
    write_trimming A_trimming_reverse.log 5 4
    write_merging B 20 18
    write_trimming B_trimming_forward.log 18 18
    write_trimming B_trimming_reverse.log 12 11
    run bash "${SCRIPT}" ids.txt out.tsv
    [ "${status}" -eq 0 ]
    run tail -n 1 out.tsv
    [ "${output}" = "$(printf 'Total\t30\t27\t27\t17\t15')" ]
}

@test "rows are sorted by sample id regardless of ids-file order" {
    printf 'B\nA\n' > ids.txt
    write_merging A 1 1
    write_trimming A_trimming_forward.log 1 1
    write_trimming A_trimming_reverse.log 1 1
    write_merging B 2 2
    write_trimming B_trimming_forward.log 2 2
    write_trimming B_trimming_reverse.log 2 2
    run bash "${SCRIPT}" ids.txt out.tsv
    [ "${status}" -eq 0 ]
    run sed -n '2p' out.tsv
    [ "${output}" = "$(printf 'A\t1\t1\t1\t1\t1')" ]
    run sed -n '3p' out.tsv
    [ "${output}" = "$(printf 'B\t2\t2\t2\t2\t2')" ]
}
