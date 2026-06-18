process build_occurrence_table {
    // [S35]: merge swarm / uchime / quality / stampa / distribution
    // into the filtered occurrence table. Representatives, stats,
    // swarms, chimera annotations, and taxonomic assignments are
    // each concatenated across the pre-cleave outputs and the
    // post-cleave ([S22]/[S36]/[S37]) outputs so cleaved clusters
    // get a row.
    //
    // Concat order matters for the taxonomic assignments: stampa_parse
    // overwrites by amplicon ID, so the file listed *last* on stdin
    // wins for an overlapping seed. Bash uses ``${TAX}{2,}`` (results2
    // first, results last) so the pre-cleave assignment wins on
    // overlap — that ordering is preserved here.
    //
    // [S46]: this is an intermediate OTU table — it is **not**
    // published; only the final <basename>_table.tsv reaches the
    // results folder.

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
    path assignments,       stageAs: 'assignments_global'
    path assignments_2,     stageAs: 'assignments_cleaved'
    path distribution
    val basename
    val sample_ids   // [S09]: comma-separated sample IDs; empty samples
                     // contribute no .distr rows but still need a zero
                     // column in the occurrence table.

    output:
    path "${basename}.OTU.filtered.cleaved.table", emit: table

    shell:
    '''
    build_filtered_contingency_table.py \
        --representatives <(cat !{representatives} !{representatives_2}) \
        --stats           <(cat !{stats} !{stats_2}) \
        --swarms          <(cat !{swarms} !{swarms_2}) \
        --chimera         <(cat !{uchime} !{uchime_2}) \
        --quality         !{quality} \
        --assignments     <(cat !{assignments_2} !{assignments}) \
        --distribution    !{distribution} \
        --samples         '!{sample_ids}' \
        > !{basename}.OTU.filtered.cleaved.table
    '''
}
