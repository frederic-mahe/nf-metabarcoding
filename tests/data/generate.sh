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

# Build a quality string of N copies of an arbitrary character. Used for
# the maximum-quality fixture ([S106]) where the fill character is '~'
# (Phred 93). Kept separate from qual_string so the SC2059-safe literal
# 'I' format above is left untouched.
qual_string_char() {
    local -r n="${1}"
    local -r char="${2}"
    local blanks
    blanks="$(printf '%*s' "${n}" '')"
    printf '%s' "${blanks// /${char}}"
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

emit_highqual_pair() {
    # [S106] a mergeable paired-end fixture whose every base carries the
    # maximum Phred+33 quality score (Q93, the '~' character). vsearch's
    # default --fastq_qmax is 41 and aborts on any higher score, so this
    # fixture only survives the fastq-reading processes when they pass
    # --fastq_qmax 93 (the range PacBio HiFi reads use). A single
    # amplicon is enough to exercise the read + merge path; because the
    # overlap of two agreeing Q93 reads has a posterior quality far above
    # Q41, the merged output also witnesses --fastq_qmaxout 93.
    local -r prefix="highqual"
    local -r r1="${prefix}_1.fastq"
    local -r r2="${prefix}_2.fastq"
    local -r amplicon="${AMPLICONS_OK[0]}"
    local fwd
    fwd="${amplicon:0:${READ_LEN}}"
    local rev_src
    rev_src="${amplicon: -${READ_LEN}}"
    local rev
    rev="$(reverse_complement "${rev_src}")"
    local qf
    qf="$(qual_string_char "${#fwd}" '~')"
    local qr
    qr="$(qual_string_char "${#rev}" '~')"
    printf '@%s\n%s\n+\n%s\n' "read_1 1:N:0:1" "${fwd}" "${qf}" > "${r1}"
    printf '@%s\n%s\n+\n%s\n' "read_1 2:N:0:1" "${rev}" "${qr}" > "${r2}"
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
    # Shadow-pipeline sample ([S04]/[S56]) — fed into the shadow
    # Part B workflow `part_B_shadow`. The sequence carries the
    # A-padded join site emitted by Part A's shadow pipeline
    # (`vsearch --fastq_join --join_padgap AAAAAAAA` by default;
    # see [S63]). A second clean record gives swarm something
    # non-trivial to cluster against.
    local -r sha_d="dddddddddddddddddddddddddddddddddddddddd"
    local -r seq_a_padded="ACGTACGTAAAAAAAAACGTACGTACGTACGTACGTACGT"
    local -r seq_d="GGGGTTTTGGGGCCCCGGGGTTTTGGGGCCCCGGGGTTTT"
    {
        printf '>%s;size=5\n%s\n' "${sha_n}" "${seq_a_padded}"
        printf '>%s;size=3\n%s\n' "${sha_d}" "${seq_d}"
    } > "${dir}/S1_notmerged.fas"
    {
        printf '%s 0.010000 40\n' "${sha_n}"
        printf '%s 0.020000 40\n' "${sha_d}"
    } > "${dir}/S1_notmerged.qual"
    {
        printf '5\t5\t%s\t5\t0\t1\t1\n' "${sha_n}"
        printf '3\t3\t%s\t3\t0\t1\t1\n' "${sha_d}"
    } > "${dir}/S1_notmerged.stats"

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

emit_compression_variants() {
    # [S03]: three byte-identical fastq pairs in different
    # compression forms — plain, gzip, bzip2. The amplicon set is
    # the same as paired_merge_ok so vsearch / cutadapt / swarm
    # have something to merge / trim / cluster.
    #
    # The three pairs share one folder so a single workflow run
    # discovers and processes them all together; distinct
    # filename prefixes (`cmp_plain_`, `cmp_gz_`, `cmp_bz2_`) keep
    # the derived sample IDs ([S12]) collision-free per [S13].
    local -r dir="compression_variants"
    mkdir -p "${dir}"

    local i
    local fwd
    local rev_src
    local rev
    : > "${dir}/cmp_plain_1.fastq"
    : > "${dir}/cmp_plain_2.fastq"
    i=0
    for A in "${AMPLICONS_OK[@]}"; do
        i=$((i + 1))
        fwd="${A:0:${READ_LEN}}"
        rev_src="${A: -${READ_LEN}}"
        rev="$(reverse_complement "${rev_src}")"
        write_record "${dir}/cmp_plain_1.fastq" "read_${i} 1:N:0:1" "${fwd}"
        write_record "${dir}/cmp_plain_2.fastq" "read_${i} 2:N:0:1" "${rev}"
    done
    # Mirror the same content into gz and bz2 siblings.
    gzip --keep --force "${dir}/cmp_plain_1.fastq"
    gzip --keep --force "${dir}/cmp_plain_2.fastq"
    mv "${dir}/cmp_plain_1.fastq.gz" "${dir}/cmp_gz_1.fastq.gz"
    mv "${dir}/cmp_plain_2.fastq.gz" "${dir}/cmp_gz_2.fastq.gz"
    bzip2 --keep --force "${dir}/cmp_plain_1.fastq"
    bzip2 --keep --force "${dir}/cmp_plain_2.fastq"
    mv "${dir}/cmp_plain_1.fastq.bz2" "${dir}/cmp_bz2_1.fastq.bz2"
    mv "${dir}/cmp_plain_2.fastq.bz2" "${dir}/cmp_bz2_2.fastq.bz2"
}


emit_hash_function_fixtures() {
    # [S65]: MD5-width (32 hex char) fixtures used to prove the .qual
    # dedup steps track params.hash_function rather than assuming the
    # 40-char SHA1 width. IDs are built from a single letter so tests
    # can assert on them readably (these stand in for md5 digests).
    local -r dir="hash_function"
    mkdir -p "${dir}"

    local -r md5_a="aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
    local -r md5_b="bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb"
    local -r md5_c="cccccccccccccccccccccccccccccccc"
    local -r seq_a="ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
    local -r seq_c="GGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTT"

    # filter_and_convert_to_fasta output shape: `>MD5;ee=<f>;length=<i>`
    # at --fasta_width 0. Two records share md5_a with different ee so
    # extract_expected_error_values must collapse them onto the lower-ee
    # row — which only works when uniq --check-chars is 32, not 40.
    {
        printf '>%s;ee=0.020000;length=40\n%s\n' "${md5_a}" "${seq_a}"
        printf '>%s;ee=0.010000;length=40\n%s\n' "${md5_a}" "${seq_a}"
        printf '>%s;ee=0.030000;length=40\n%s\n' "${md5_c}" "${seq_c}"
    } > "${dir}/md5_filtered.fas"

    # Per-sample .qual (extract_ee.awk format) with 32-char names, for
    # build_expected_error_file. S2 carries the lower ee for md5_a so
    # the merge must keep S2's row.
    {
        printf '%s 0.010000 40\n' "${md5_a}"
        printf '%s 0.020000 40\n' "${md5_b}"
    } > "${dir}/S1_md5.qual"
    {
        printf '%s 0.005000 40\n' "${md5_a}"
        printf '%s 0.030000 40\n' "${md5_c}"
    } > "${dir}/S2_md5.qual"
}

emit_large_cleaved_fixture() {
    # [S81]: a cleaved-representatives fasta large enough that the
    # `sed ... | sort -n | head -n 1` minsize computation in
    # chimera_detection_post_cleave overflows the OS pipe buffer
    # (~64 KB of sorted size numbers) and SIGPIPEs `sort` under
    # `set -euo pipefail`. 10k records with 6-digit sizes (~70 KB of
    # sorted output) reliably trigger it; the SIGPIPE-safe `sed -n 1p`
    # must survive. Sequences are identical so the downstream uchime
    # pass stays fast.
    local -r dir="part_b"
    mkdir -p "${dir}"
    awk 'BEGIN {
        seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
        for (i = 1; i <= 10000; i++)
            printf ">cleaved_%d;size=%d\n%s\n", i, 100000 + i, seq
    }' > "${dir}/big_cleaved.fas2"
}

emit_recluster_fixture() {
    # [S102]-[S105]/D20: per-sample Part B inputs for the optional
    # post-mumu re-clustering pass. Three 40-nt variants of one seed,
    # each two substitutions from the seed (so swarm d=1 keeps them as
    # three separate OTUs and fastidious grafting — mass >= boundary 3 —
    # does not bridge them) and each dominant in a distinct sample (so
    # they never co-occur and mumu keeps all three). They are 90-95%
    # identical, so `vsearch --cluster_size --id 0.9` folds all three
    # into one centroid: the OFF run emits three rows, the ON run one.
    local -r dir="recluster"
    mkdir -p "${dir}"

    local -r id0="1111111111111111111111111111111111111111"
    local -r id1="2222222222222222222222222222222222222222"
    local -r id2="3333333333333333333333333333333333333333"

    # v0 = seed; v1 = seed with a 3' CG->AA; v2 = seed with a 5' GA->TT.
    local -r v0="GATCAGTCAGTCAGGTCAGTGCATGCATGCATTAGCATCG"
    local -r v1="GATCAGTCAGTCAGGTCAGTGCATGCATGCATTAGCATAA"
    local -r v2="TTTCAGTCAGTCAGGTCAGTGCATGCATGCATTAGCATCG"

    # Per-sample .fas (vsearch "SHA1;size=N" header, --fasta_width 0).
    printf '>%s;size=20\n%s\n' "${id0}" "${v0}" > "${dir}/rcA.fas"
    printf '>%s;size=10\n%s\n' "${id1}" "${v1}" > "${dir}/rcB.fas"
    printf '>%s;size=8\n%s\n'  "${id2}" "${v2}" > "${dir}/rcC.fas"

    # Per-sample .qual (extract_ee.awk "<SHA1> <ee> <length>"). ee/length
    # = 0.0001 <= max_ee 0.0002 so every variant clears the [S35] filter.
    printf '%s 0.004000 40\n' "${id0}" > "${dir}/rcA.qual"
    printf '%s 0.004000 40\n' "${id1}" > "${dir}/rcB.qual"
    printf '%s 0.004000 40\n' "${id2}" > "${dir}/rcC.qual"

    # Per-sample .stats (swarm --statistics-file rows). One single-
    # amplicon local cluster per sample.
    printf '1\t20\t%s\t20\t1\t0\t0\n' "${id0}" > "${dir}/rcA.stats"
    printf '1\t10\t%s\t10\t1\t0\t0\n' "${id1}" > "${dir}/rcB.stats"
    printf '1\t8\t%s\t8\t1\t0\t0\n'  "${id2}" > "${dir}/rcC.stats"
}


emit_paired "paired_merge_ok"   "${AMPLICONS_OK[@]}"
emit_paired "paired_merge_fail" "${AMPLICONS_LONG[@]}"
emit_highqual_pair
emit_single
emit_uncompressed
emit_empty
emit_unpaired_only
emit_miseq_pair
emit_reserved_keyword
emit_duplicate_sample_ids
emit_part_b_fixtures
emit_e2e_part_b_fixture
emit_compression_variants
emit_hash_function_fixtures
emit_large_cleaved_fixture
emit_recluster_fixture

echo "Wrote fixtures to ${DATA_DIR}"
