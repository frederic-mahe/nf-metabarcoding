include { normalize_path } from '../functions.nf'


process merge_substring_otus {
    // [S39]: wraps merge_sub_superstring_OTUs_with_larger_OTUs.py
    // + the legacy bash `sort_occurrence_table` step
    // + the read-count invariant check. Pupil OTUs are merged into
    // their masters (samples summed, spread/total/cloud updated),
    // the resulting table is sorted by the OTU column, and we
    // assert that the global read count is conserved.
    //
    // [S45]: cat the upstream vsearch search.log with the merge
    // step's stderr to produce the combined
    // <basename>_superstring_clustering.log. The merged OTU table
    // itself is **not** published ([S46]).
    publishDir path: { normalize_path(params.results_folder) }, mode: params.publish_mode,
        pattern: "*_superstring_clustering.log",
        enabled: params.results_folder != null

    input:
    path otu_table
    path matches
    path search_log
    val basename

    output:
    path "${otu_table.baseName}.nosubstringOTUs.table", emit: table
    path "${basename}_superstring_clustering.log",      emit: log

    shell:
    '''
    #!/bin/bash
    set -euo pipefail

    out="!{otu_table.baseName}.nosubstringOTUs.table"
    log="!{basename}_superstring_clustering.log"
    tmp_table="$(mktemp)"
    merge_stderr="$(mktemp)"
    trap 'rm -f "${tmp_table}" "${merge_stderr}"' EXIT

    {
        echo "=== search_for_terminal_gaps (vsearch --cluster_smallmem) ==="
        cat !{search_log}
        echo
        echo "=== merge_substring_otus + invariant check ==="
    } > "${log}"

    {
        merge_substring_otus.py \
            -t !{otu_table} \
            -m !{matches} \
            -o "${tmp_table}"

        (head -n 1 "${tmp_table}"
         tail -n +2 "${tmp_table}" | sort -k1,1n) > "${out}"

        before="$(awk 'NR > 1 {t += $2} END {print t + 0}' !{otu_table})"
        after="$(awk 'NR > 1 {t += $2} END {print t + 0}' "${out}")"
        if (( before != after )) ; then
            echo "merge_substring_otus: read count changed (${before} -> ${after})" >&2
            exit 1
        fi
        echo "merge_substring_otus: read count conserved (${before})" >&2
    } 2> "${merge_stderr}"
    cat "${merge_stderr}" >> "${log}"
    '''
}
