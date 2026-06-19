#!/usr/bin/env bash
# [S76] Generate the bundled demo dataset used by `-profile demo`.
#
# Unlike tests/data/generate.sh (whose fastq outputs are .gitignored and
# regenerated per test run), this script's outputs are **committed** as
# plain text so `nextflow run main.nf -profile demo` works in a fresh
# checkout with no setup. Re-run it only to regenerate the committed
# files (deterministic; safe to re-run):
#
#     bash assets/demo/make_demo.sh
#
# The reads are synthetic and tiny — enough to exercise the whole
# pipeline (merge -> trim -> dereplicate -> cluster -> occurrence table
# -> taxonomy) so a site can confirm its container / cluster setup, not
# to mean anything biologically.
set -euo pipefail

OUT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
readonly OUT_DIR
cd "${OUT_DIR}"

# Primers match the demo-profile defaults (the IUPAC ambiguities of the
# real defaults collapsed to a single base, as in tests/data/generate.sh:
# forward CCAGCASCYGCGGTAATTCC, reverse ACTTTCGTTCTTGATYRA).
readonly FORWARD_PRIMER="CCAGCACCCGCGGTAATTCC"
readonly REVERSE_PRIMER_RC="TTGATCAAGAACGAAAGT"

# Two distinct 60 nt amplicon bodies (reused from generate.sh). Each
# becomes one OTU.
readonly BODY1="GATCAGATACCGTCGTAGTCTTAACCATAAACTATGCCGACTAGGGATCGGGCGATGTTA"
readonly BODY2="CGTACGATCCAGCATTGGATCCATAATAGTCATCAAGTCAGGGATCGTCCTGAATCGGAT"

readonly READ_LEN=70
# Copies per amplicon. >2 so each OTU clears the per-cluster minimum
# ([S17]) and lands a row in the occurrence table.
readonly COPIES=5

reverse_complement() {
    local -r nucleotides="acgturykmbdhvswACGTURYKMBDHVSW"
    local -r complements="tgcaayrmkvhdbswTGCAAYRMKVHDBSW"
    tr "${nucleotides}" "${complements}" <<< "${1}" | rev
}

qual_string() {
    local -r n="${1}"
    printf 'I%.0s' $(seq 1 "${n}")
}

write_record() {
    local -r out="${1}"
    local -r name="${2}"
    local -r seq="${3}"
    local q
    q="$(qual_string "${#seq}")"
    printf '@%s\n%s\n+\n%s\n' "${name}" "${seq}" "${q}" >> "${out}"
}

readonly R1="demo_1.fastq"
readonly R2="demo_2.fastq"
: > "${R1}"
: > "${R2}"

emit_pairs() {
    local -r body="${1}"
    local -r tag="${2}"
    local amplicon
    amplicon="${FORWARD_PRIMER}${body}${REVERSE_PRIMER_RC}"
    local fwd
    fwd="${amplicon:0:${READ_LEN}}"
    local rev
    rev="$(reverse_complement "${amplicon: -${READ_LEN}}")"
    local i
    for (( i = 1; i <= COPIES; i++ )); do
        write_record "${R1}" "${tag}_${i} 1:N:0:1" "${fwd}"
        write_record "${R2}" "${tag}_${i} 2:N:0:1" "${rev}"
    done
}

emit_pairs "${BODY1}" "otu1"
emit_pairs "${BODY2}" "otu2"

# Stampa-formatted reference ([S47]: `>id <space> |-separated lineage`).
# The two reference sequences are the trimmed amplicon bodies, so the
# representatives match at 100 % identity and Part C assigns a lineage.
cat > reference.fasta <<REF
>demo_ref_1 Bacteria|Firmicutes|Bacilli|Lactobacillales|Lactobacillaceae|Lactobacillus
${BODY1}
>demo_ref_2 Bacteria|Proteobacteria|Gammaproteobacteria|Enterobacterales|Enterobacteriaceae|Escherichia
${BODY2}
REF

echo "wrote ${R1}, ${R2}, reference.fasta to ${OUT_DIR}"
