include { log_dir }     from '../functions.nf'
include { publish_dir } from '../functions.nf'


process recluster_merge {
    // [S104]/[S105]/D20: folds each re-clustered member OTU onto its
    // vsearch centroid via bin/recluster_otu_table.py (samples + total
    // summed, spread recomputed, metadata from the centroid, cloud left
    // at "NA", OTUs renumbered 1..N), then asserts the global read count
    // is conserved (like [S39]).
    //
    // [S105]/D-a: when --recluster_id is set the reclustered table
    // *replaces* the post-mumu table as Part B's deliverable, so this
    // process publishes the final <basename>_table.tsv into
    // occurrence_table/ ([S59]/[S71]) — the fine pre-recluster table
    // (rebuild_post_mumu_table) is left unpublished. The post-mumu table
    // is staged under a fixed name so the <basename>_table.tsv output
    // never collides with the identically-named input.
    //
    // [S105]: the combined vsearch-search + merge stderr is published as
    // a conditional 7th step log <basename>_reclustering.log under
    // logs/part_b/ ([S45] family, D16).
    publishDir path: { publish_dir('occurrence_table') }, mode: params.publish_mode,
        pattern: "*_table.tsv"
    publishDir path: { log_dir('part_b') }, mode: params.publish_mode,
        pattern: "*_reclustering.log"

    input:
    path post_mumu_table, stageAs: 'post_mumu_table.tsv'
    path matches
    path search_log
    val basename

    output:
    path "${basename}_table.tsv",         emit: table
    path "${basename}_reclustering.log",  emit: log

    shell:
    '''
    #!/bin/bash
    set -euo pipefail

    out="!{basename}_table.tsv"
    log="!{basename}_reclustering.log"
    merge_stderr="$(mktemp)"
    trap 'rm -f "${merge_stderr}"' EXIT

    {
        echo "=== recluster_search (vsearch --cluster_size) ==="
        cat !{search_log}
        echo
        echo "=== recluster_otu_table + invariant check ==="
    } > "${log}"

    {
        recluster_otu_table.py \
            -t post_mumu_table.tsv \
            -m !{matches} \
            -o "${out}"

        before="$(awk 'NR > 1 {t += $2} END {print t + 0}' post_mumu_table.tsv)"
        after="$(awk 'NR > 1 {t += $2} END {print t + 0}' "${out}")"
        if (( before != after )) ; then
            echo "recluster_merge: read count changed (${before} -> ${after})" >&2
            exit 1
        fi
        echo "recluster_merge: read count conserved (${before})" >&2
    } 2> "${merge_stderr}"
    cat "${merge_stderr}" >> "${log}"
    '''

    stub:
    """
    touch ${basename}_table.tsv ${basename}_reclustering.log
    """
}
