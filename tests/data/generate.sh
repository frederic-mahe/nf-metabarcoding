#!/usr/bin/env bash
# Generate minimal fastq fixtures for the test suite.
# Deterministic; safe to re-run. Outputs land in tests/data/.

set -euo pipefail

DATA_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${DATA_DIR}"

# Primers must match the defaults in main.nf.
# Forward primer in main.nf is CCAGCASCYGCGGTAATTCC; we collapse the
# IUPAC ambiguities (S=C/G -> C, Y=C/T -> C) for the synthetic reads.
# Reverse primer is ACTTTCGTTCTTGATYRA (Y=C/T -> C, R=A/G -> A), so its
# reverse complement embedded at the 3' end of the amplicon is:
FORWARD_PRIMER="CCAGCACCCGCGGTAATTCC"
REVERSE_PRIMER_RC="TTGATCAAGAACGAAAGT"

# Reverse complement a DNA string.
rc() {
    tr 'ACGTacgt' 'TGCAtgca' <<< "${1}" | rev
}

# Build a quality string of N copies of 'I' (Phred 40).
qual_string() {
    local n="${1}"
    printf 'I%.0s' $(seq 1 "${n}")
}

# Write one fastq record to a file.
write_record() {
    local out="${1}" name="${2}" seq="${3}"
    local q; q="$(qual_string "${#seq}")"
    printf '@%s\n%s\n+\n%s\n' "${name}" "${seq}" "${q}" >> "${out}"
}

# Synthetic amplicon bodies (60 nt each). Combined with primers,
# each amplicon is 98 nt long (20 + 60 + 18). Bodies are non-periodic
# so paired reads have a single unambiguous overlap.
BODY1="GATCAGATACCGTCGTAGTCTTAACCATAAACTATGCCGACTAGGGATCGGGCGATGTTA"
BODY2="CGTACGATCCAGCATTGGATCCATAATAGTCATCAAGTCAGGGATCGTCCTGAATCGGAT"
BODY3="TGAATCCAGGTAACGGATCCGGATCTAGCAGTACATCGATTCAGGTACATCGACCATGAT"

AMPLICONS_OK=(
    "${FORWARD_PRIMER}${BODY1}${REVERSE_PRIMER_RC}"
    "${FORWARD_PRIMER}${BODY2}${REVERSE_PRIMER_RC}"
    "${FORWARD_PRIMER}${BODY3}${REVERSE_PRIMER_RC}"
)

# For paired_merge_fail: longer body so R1 and R2 cannot overlap with
# a 70 nt read length (insert = 220 nt > 2 * read_length).
LONG_BODY="$(printf 'ACGT%.0s' $(seq 1 50))"
AMPLICONS_LONG=(
    "${FORWARD_PRIMER}${LONG_BODY}${REVERSE_PRIMER_RC}"
)

READ_LEN=70

emit_paired() {
    local prefix="${1}" amps=("${@:2}")
    local r1="${prefix}_1.fastq" r2="${prefix}_2.fastq"
    : > "${r1}"
    : > "${r2}"
    local i=0
    for A in "${amps[@]}"; do
        i=$((i + 1))
        local fwd="${A:0:${READ_LEN}}"
        local rev_src="${A: -${READ_LEN}}"
        local rev; rev="$(rc "${rev_src}")"
        write_record "${r1}" "read_${i} 1:N:0:1" "${fwd}"
        write_record "${r2}" "read_${i} 2:N:0:1" "${rev}"
    done
    gzip --force "${r1}" "${r2}"
}

emit_single() {
    local out="single_end.fastq"
    : > "${out}"
    local i=0
    for A in "${AMPLICONS_OK[@]}"; do
        i=$((i + 1))
        write_record "${out}" "read_${i}" "${A:0:${READ_LEN}}"
    done
    gzip --force "${out}"
}

emit_uncompressed() {
    # Same content as paired_merge_ok but left uncompressed.
    local r1="uncompressed_1.fastq" r2="uncompressed_2.fastq"
    : > "${r1}"
    : > "${r2}"
    local i=0
    for A in "${AMPLICONS_OK[@]}"; do
        i=$((i + 1))
        local fwd="${A:0:${READ_LEN}}"
        local rev_src="${A: -${READ_LEN}}"
        local rev; rev="$(rc "${rev_src}")"
        write_record "${r1}" "read_${i} 1:N:0:1" "${fwd}"
        write_record "${r2}" "read_${i} 2:N:0:1" "${rev}"
    done
}

emit_empty() {
    # A valid empty gzipped fastq pair.
    : | gzip > "empty_1.fastq.gz"
    : | gzip > "empty_2.fastq.gz"
}

emit_paired "paired_merge_ok"   "${AMPLICONS_OK[@]}"
emit_paired "paired_merge_fail" "${AMPLICONS_LONG[@]}"
emit_single
emit_uncompressed
emit_empty

echo "Wrote fixtures to ${DATA_DIR}"
