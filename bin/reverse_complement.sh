#!/usr/bin/env bash
#
# reverse_complement.sh
#
# Print the reverse complement of an IUPAC nucleotide string. Handles
# the full IUPAC alphabet in both cases. N and I are their own
# complements and pass through unchanged.
#
# Usage:
#   reverse_complement.sh ACGTRY    # argument form
#   echo ACGTRY | reverse_complement.sh   # stdin form

set -euo pipefail

readonly NUCLEOTIDES="acgturykmbdhvswACGTURYKMBDHVSW"
readonly COMPLEMENTS="tgcaayrmkvhdbswTGCAAYRMKVHDBSW"

if [[ "$#" -gt 0 ]]; then
    tr "${NUCLEOTIDES}" "${COMPLEMENTS}" <<< "${1}" | rev
else
    tr "${NUCLEOTIDES}" "${COMPLEMENTS}" | rev
fi
