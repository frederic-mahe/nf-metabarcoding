#!/usr/bin/env bash
# Generate minimal fastq fixtures for the test suite.
# Deterministic; safe to re-run. Outputs land in tests/data/.

set -euo pipefail

DATA_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
readonly DATA_DIR
cd "${DATA_DIR}"

# Primers must match the defaults in main.nf.
# Forward primer in main.nf is CCAGCASCYGCGGTAATTCC; we collapse the
# IUPAC ambiguities (S=C/G -> C, Y=C/T -> C) for the synthetic reads.
# Reverse primer is ACTTTCGTTCTTGATYRA (Y=C/T -> C, R=A/G -> A), so its
# reverse complement embedded at the 3' end of the amplicon is:
readonly FORWARD_PRIMER="CCAGCACCCGCGGTAATTCC"
readonly REVERSE_PRIMER_RC="TTGATCAAGAACGAAAGT"

reverse_complement() {
    # reverse-complement a DNA/RNA IUPAC string
    # Note: N and I are their own complements, no need to include them
    local -r nucleotides="acgturykmbdhvswACGTURYKMBDHVSW"
    local -r complements="tgcaayrmkvhdbswTGCAAYRMKVHDBSW"

    tr "${nucleotides}" "${complements}" <<< "${1}" | rev
}

# Build a quality string of N copies of 'I' (Phred 40).
qual_string() {
    local -r n="${1}"
    printf 'I%.0s' $(seq 1 "${n}")
}

# Write one fastq record to a file.
write_record() {
    local -r out="${1}"
    local -r name="${2}"
    local -r seq="${3}"
    local q
    q="$(qual_string "${#seq}")"
    printf '@%s\n%s\n+\n%s\n' "${name}" "${seq}" "${q}" >> "${out}"
}

# Synthetic amplicon bodies (60 nt each). Combined with primers,
# each amplicon is 98 nt long (20 + 60 + 18). Bodies are non-periodic
# so paired reads have a single unambiguous overlap.
readonly BODY1="GATCAGATACCGTCGTAGTCTTAACCATAAACTATGCCGACTAGGGATCGGGCGATGTTA"
readonly BODY2="CGTACGATCCAGCATTGGATCCATAATAGTCATCAAGTCAGGGATCGTCCTGAATCGGAT"
readonly BODY3="TGAATCCAGGTAACGGATCCGGATCTAGCAGTACATCGATTCAGGTACATCGACCATGAT"

AMPLICONS_OK=(
    "${FORWARD_PRIMER}${BODY1}${REVERSE_PRIMER_RC}"
    "${FORWARD_PRIMER}${BODY2}${REVERSE_PRIMER_RC}"
    "${FORWARD_PRIMER}${BODY3}${REVERSE_PRIMER_RC}"
)
readonly -a AMPLICONS_OK

# For paired_merge_fail: longer body so R1 and R2 cannot overlap with
# a 70 nt read length (insert = 220 nt > 2 * read_length).
LONG_BODY="$(printf 'ACGT%.0s' $(seq 1 50))"
readonly LONG_BODY

AMPLICONS_LONG=(
    "${FORWARD_PRIMER}${LONG_BODY}${REVERSE_PRIMER_RC}"
)
readonly -a AMPLICONS_LONG

readonly READ_LEN=70

emit_paired() {
    local -r prefix="${1}"
    local -ra amps=("${@:2}")
    local -r r1="${prefix}_1.fastq"
    local -r r2="${prefix}_2.fastq"
    : > "${r1}"
    : > "${r2}"
    local i=0
    for A in "${amps[@]}"; do
        i=$((i + 1))
        local fwd
        local rev_src
        local rev
        fwd="${A:0:${READ_LEN}}"
        rev_src="${A: -${READ_LEN}}"
        rev="$(reverse_complement "${rev_src}")"
        write_record "${r1}" "read_${i} 1:N:0:1" "${fwd}"
        write_record "${r2}" "read_${i} 2:N:0:1" "${rev}"
    done
    gzip --force "${r1}" "${r2}"
}

emit_single() {
    local -r out="single_end.fastq"
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
    local -r r1="uncompressed_1.fastq"
    local -r r2="uncompressed_2.fastq"
    : > "${r1}"
    : > "${r2}"
    local i=0
    for A in "${AMPLICONS_OK[@]}"; do
        i=$((i + 1))
        local fwd
        local rev_src
        local rev
        fwd="${A:0:${READ_LEN}}"
        rev_src="${A: -${READ_LEN}}"
        rev="$(reverse_complement "${rev_src}")"
        write_record "${r1}" "read_${i} 1:N:0:1" "${fwd}"
        write_record "${r2}" "read_${i} 2:N:0:1" "${rev}"
    done
}

emit_empty() {
    # A valid empty gzipped fastq pair.
    : | gzip > "empty_1.fastq.gz"
    : | gzip > "empty_2.fastq.gz"
}

emit_unpaired_only() {
    # An isolated directory containing a single full-length single-end
    # fastq file: name does not match any paired-end pattern, so the
    # workflow must route it through the unpaired branch ([S21]).
    local -r dir="unpaired_only"
    mkdir -p "${dir}"
    local -r out="${dir}/unpaired_sample.fastq"
    : > "${out}"
    local i=0
    for A in "${AMPLICONS_OK[@]}"; do
        i=$((i + 1))
        write_record "${out}" "read_${i}" "${A}"
    done
    gzip --force "${out}"
}

emit_miseq_pair() {
    # A MiSeq-style paired-end fixture exercising canonical pattern
    # row 1 (`_L00[1-9]_R1_00[1-9]\.<ext>`). The sample ID derived
    # from the R1 basename is `SampleX_S1`.
    local -r dir="miseq"
    mkdir -p "${dir}"
    local -r r1="${dir}/SampleX_S1_L001_R1_001.fastq"
    local -r r2="${dir}/SampleX_S1_L001_R2_001.fastq"
    : > "${r1}"
    : > "${r2}"
    local i=0
    for A in "${AMPLICONS_OK[@]}"; do
        i=$((i + 1))
        local fwd
        local rev_src
        local rev
        fwd="${A:0:${READ_LEN}}"
        rev_src="${A: -${READ_LEN}}"
        rev="$(reverse_complement "${rev_src}")"
        write_record "${r1}" "read_${i} 1:N:0:1" "${fwd}"
        write_record "${r2}" "read_${i} 2:N:0:1" "${rev}"
    done
    gzip --force "${r1}" "${r2}"
}

emit_paired "paired_merge_ok"   "${AMPLICONS_OK[@]}"
emit_paired "paired_merge_fail" "${AMPLICONS_LONG[@]}"
emit_single
emit_uncompressed
emit_empty
emit_unpaired_only
emit_miseq_pair

echo "Wrote fixtures to ${DATA_DIR}"
