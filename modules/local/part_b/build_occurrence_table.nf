process build_occurrence_table {
    // [S35]: merge swarm / uchime / quality / distribution into the
    // filtered occurrence table. Representatives, stats, swarms, and
    // chimera annotations are each concatenated across the pre-cleave
    // outputs and the post-cleave ([S22]/[S37]) outputs so cleaved
    // clusters get a row.
    //
    // [S127]: no taxonomic assignments exist yet at this point, so
    // --assignments is not passed and the script writes the
    // placeholder identity / taxonomy / references (0.0 / NA / NA).
    // Part C's update_occurrence_table overwrites them later.
    //
    // [S46]: this is an intermediate OTU table — it is **not**
    // published; only the final <basename>_table.tsv reaches the
    // results folder.
    //
    // [S135]: the script's stderr (progress lines and the final
    // FILTER: summary) is kept in <basename>_occurrence_filter.log,
    // emitted as filter_log and not published (the [S45] logs/part_b/
    // file list stays as is). Part B reads it to warn about an empty
    // table. On failure the log is copied to stderr, so the error
    // still shows in Nextflow's report.

    input:
    path representatives,   stageAs: 'reps_global.fas'
    path representatives_2, stageAs: 'reps_cleaved.fas2'
    path stats,             stageAs: 'stats_global'
    path stats_2,           stageAs: 'stats_cleaved'
    path swarms,            stageAs: 'swarms_global'
    path swarms_2,          stageAs: 'swarms_cleaved'
    path uchime,            stageAs: 'uchime_global'
    path uchime_2,          stageAs: 'uchime_cleaved'
    path quality
    path distribution
    val basename
    val sample_ids   // [S09]: comma-separated sample IDs; empty samples
                     // contribute no .distr rows but still need a zero
                     // column in the occurrence table.

    output:
    path "${basename}.OTU.filtered.cleaved.table", emit: table
    path "${basename}_occurrence_filter.log",      emit: filter_log

    shell:
    '''
    readonly FILTER_LOG="!{basename}_occurrence_filter.log"

    build_filtered_contingency_table.py \
        --representatives <(cat !{representatives} !{representatives_2}) \
        --stats           <(cat !{stats} !{stats_2}) \
        --swarms          <(cat !{swarms} !{swarms_2}) \
        --chimera         <(cat !{uchime} !{uchime_2}) \
        --quality         !{quality} \
        --distribution    !{distribution} \
        --samples         '!{sample_ids}' \
        --max-ee          !{params.max_ee} \
        --min-abundance   !{params.min_abundance} \
        --min-spread      !{params.min_spread} \
        > !{basename}.OTU.filtered.cleaved.table \
        2> "${FILTER_LOG}" || {
        status=$?
        cat "${FILTER_LOG}" >&2
        exit "${status}"
    }
    '''

    stub:
    """
    touch ${basename}.OTU.filtered.cleaved.table ${basename}_occurrence_filter.log
    """
}
