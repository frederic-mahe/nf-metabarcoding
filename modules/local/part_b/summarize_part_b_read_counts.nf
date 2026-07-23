include { log_dir } from '../functions.nf'


process summarize_part_b_read_counts {
    // [S107]: per-sample read/cluster tracking summary for Part B — the
    // counterpart of Part A's summarize_read_counts ([S86]). Part B is
    // pooled (every sample is merged at global_dereplication, [S31]), so
    // its step logs ([S45]) carry no per-sample dimension; the summary is
    // reconstructed by bin/build_part_b_read_counts.py from the .distr
    // ([S29]) and the per-sample columns of the intermediate OTU tables.
    // Published to logs/part_b/ ([D16]) beside the step logs, one table
    // per run (a project-wide Part B summary).
    //
    // The `final_table` input is the emitted table ([S46]/[S105]): the
    // post-mumu table when re-clustering is off (then it is the same file
    // as `mumu_table`, staged a second time — harmless) and the
    // reclustered table when on. It is only read when `recluster` is
    // true, in which case a `clusters_recluster` column is appended.
    publishDir path: { log_dir('part_b') }, mode: params.publish_mode, pattern: "*.tsv"

    input:
    path distr,          stageAs: 'distribution.distr'
    path filtered_table, stageAs: 'filtered.table'
    path merged_table,   stageAs: 'merged.table'
    path mumu_table,     stageAs: 'mumu.table'
    path final_table,    stageAs: 'final.table'
    val  sample_ids      // [S09]: comma-separated authoritative sample IDs
    val  recluster       // boolean: append the clusters_recluster column
    val  basename

    output:
    path "${basename}_read_counts.tsv", emit: table

    shell:
    '''
    set -euo pipefail

    recluster_arg=""
    if [ "!{recluster}" = "true" ] ; then
        recluster_arg="--recluster final.table"
    fi

    build_part_b_read_counts.py \
        --samples  '!{sample_ids}' \
        --distr    distribution.distr \
        --filtered filtered.table \
        --merged   merged.table \
        --mumu     mumu.table \
        ${recluster_arg} \
        -o !{basename}_read_counts.tsv
    '''

    stub:
    """
    touch ${basename}_read_counts.tsv
    """
}
