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

emit_e2e_part_b_fixture() {
    # Isolated single-pair folders for the Part A → Part B end-to-end
    # tests. Sharing tests/data/ with the rest of the fixtures would
    # inflate N (every paired-end fastq there would be discovered as a
    # separate sample), so the e2e tests need their own clean folders.
    local -r ok="e2e_part_b/ok"
    local -r fail="e2e_part_b/fail"
    mkdir -p "${ok}" "${fail}"
    cp paired_merge_ok_1.fastq.gz   "${ok}/sampleA_1.fastq.gz"
    cp paired_merge_ok_2.fastq.gz   "${ok}/sampleA_2.fastq.gz"
    cp paired_merge_fail_1.fastq.gz "${fail}/sampleB_1.fastq.gz"
    cp paired_merge_fail_2.fastq.gz "${fail}/sampleB_2.fastq.gz"
}

emit_reserved_keyword() {
    # [S23] — a paired-end fixture whose R1 name resolves to the
    # reserved sample ID `X_notmerged`. Isolated in its own folder so
    # the validation can be exercised end-to-end without contaminating
    # the main `tests/data/` listing.
    local -r dir="reserved_keyword"
    mkdir -p "${dir}"
    cp paired_merge_ok_1.fastq.gz "${dir}/X_notmerged_1.fastq.gz"
    cp paired_merge_ok_2.fastq.gz "${dir}/X_notmerged_2.fastq.gz"
}

emit_part_b_fixtures() {
    # Hand-crafted per-sample Part A outputs for Part B's process
    # tests. Two samples (S1, S2) plus a shadow-pipeline artefact
    # (S1_notmerged.fas) that must be filtered out of the fasta
    # channel ([S27]).
    #
    # Sequence A is shared between S1 (size=3) and S2 (size=1) so
    # global_dereplication has something to collapse; sequences B
    # (S1) and C (S2) are unique to their sample.
    #
    # SHA1-shaped IDs are constructed from a single letter so tests
    # can assert on them readably.
    local -r dir="part_b"
    mkdir -p "${dir}"

    local -r sha_a="aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
    local -r sha_b="bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb"
    local -r sha_c="cccccccccccccccccccccccccccccccccccccccc"
    local -r sha_n="9999999999999999999999999999999999999999"

    local -r seq_a="ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
    local -r seq_b="TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCC"
    local -r seq_c="GGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTT"
    local -r seq_n="NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"

    # Per-sample .fas — vsearch-style "SHA1;size=N" header, no
    # trailing semicolon, --fasta_width 0 (one-line sequence).
    {
        printf '>%s;size=3\n%s\n' "${sha_a}" "${seq_a}"
        printf '>%s;size=2\n%s\n' "${sha_b}" "${seq_b}"
    } > "${dir}/S1.fas"
    {
        printf '>%s;size=1\n%s\n' "${sha_a}" "${seq_a}"
        printf '>%s;size=4\n%s\n' "${sha_c}" "${seq_c}"
    } > "${dir}/S2.fas"
    # Shadow-pipeline output — must be excluded from the Part B
    # channel ([S27]). Content is intentionally bogus so a leak
    # would surface in dedup / distribution assertions.
    {
        printf '>%s;size=99\n%s\n' "${sha_n}" "${seq_n}"
    } > "${dir}/S1_notmerged.fas"

    # Per-sample .qual — extract_ee.awk format
    # "<SHA1> <ee> <length>", sorted by length / SHA1 / ee. S2 carries
    # a lower ee for sequence A so the merge picks S2's row.
    {
        printf '%s 0.010000 40\n' "${sha_a}"
        printf '%s 0.020000 40\n' "${sha_b}"
    } > "${dir}/S1.qual"
    {
        printf '%s 0.005000 40\n' "${sha_a}"
        printf '%s 0.030000 40\n' "${sha_c}"
    } > "${dir}/S2.qual"

    # Per-sample .stats — swarm --statistics-file rows
    # (unique amplicons, total abundance, seed, seed-abundance,
    # singletons, max-generation, max-radius). Two rows in S1, one
    # in S2.
    {
        printf '3\t5\t%s\t3\t0\t1\t1\n' "${sha_a}"
        printf '2\t4\t%s\t2\t0\t1\t1\n' "${sha_b}"
    } > "${dir}/S1.stats"
    {
        printf '4\t6\t%s\t4\t0\t1\t1\n' "${sha_c}"
    } > "${dir}/S2.stats"

    # E — an empty sample. [S09]/[S27]: the empty .fas must travel
    # through to the occurrence table (downstream processes must
    # tolerate zero-record inputs). The matching .qual and .stats
    # are also empty (a sample with no surviving reads).
    : > "${dir}/E.fas"
    : > "${dir}/E.qual"
    : > "${dir}/E.stats"
}

emit_duplicate_sample_ids() {
    # [S13]/[S14]: two folders that both resolve to the same sample
    # ID `A` (one via canonical pattern row 5, the other via row 6),
    # so the workflow can verify the abort-on-duplicate behaviour.
    local -r root="duplicate_sample_ids"
    local -r a="${root}/run1"
    local -r b="${root}/run2"
    mkdir -p "${a}" "${b}"
    cp paired_merge_ok_1.fastq.gz "${a}/A_1.fastq.gz"
    cp paired_merge_ok_2.fastq.gz "${a}/A_2.fastq.gz"
    cp paired_merge_ok_1.fastq.gz "${b}/A_R1.fastq.gz"
    cp paired_merge_ok_2.fastq.gz "${b}/A_R2.fastq.gz"
}

emit_n_containing_fasta() {
    # [S04] — minimal fixture for `mask_ns_for_swarm`: one sequence
    # with 8 internal Ns (mirrors the `vsearch --fastq_join` padding),
    # one clean sequence. Sequence lines must be rewritten N→A;
    # header lines (including the literal 'N' in '>seq1_with_Ns') must
    # be left alone.
    cat > n_containing.fas <<'FASTA'
>seq1_with_Ns
ACGTACGTNNNNNNNNACGTACGT
>seq2_no_Ns
ACGTACGT
FASTA
}

emit_u_containing_fastq() {
    # [S52] — fixture for `filter_and_convert_to_fasta`'s U→T
    # normalisation. Two records (length 40, above the
    # MIN_LENGTH=32 cut-off): the first uses uppercase `U`, the
    # second uses lowercase `u`, so the test can pin both
    # branches of the case-preserving strip. Quality string is
    # Phred-40 ('I' x 40) — the literal 'U' / 'u' must NOT appear
    # on the sequence line of the resulting fasta.
    local qual
    qual="$(qual_string 40)"
    {
        printf '@u_upper\nACGUACGUACGUACGUACGUACGUACGUACGUACGUACGU\n+\n%s\n' "${qual}"
        printf '@u_lower\nacguacguacguacguacguacguacguacguacguacgu\n+\n%s\n' "${qual}"
    } > u_containing.fastq
}

emit_paired "paired_merge_ok"   "${AMPLICONS_OK[@]}"
emit_paired "paired_merge_fail" "${AMPLICONS_LONG[@]}"
emit_single
emit_uncompressed
emit_empty
emit_unpaired_only
emit_miseq_pair
emit_reserved_keyword
emit_duplicate_sample_ids
emit_part_b_fixtures
emit_e2e_part_b_fixture
emit_n_containing_fasta
emit_u_containing_fastq

echo "Wrote fixtures to ${DATA_DIR}"
