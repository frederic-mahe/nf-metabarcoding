process list_all_cluster_seeds_of_size_greater_than_2 {
    // [S30]: concatenate every per-sample <sampleId>.stats (already
    // filtered to clusters > 2 reads by Part A's list_local_clusters)
    // into a single project-wide file. Each row is prefixed with the
    // sample ID derived from the .stats filename.
    //
    // [S59]: this aggregated per-sample-OTUs stats file is an
    // internal intermediate consumed by cleaving — not published.

    input:
    path stats_files
    val basename

    output:
    path "${basename}_per_sample_OTUs.stats"

    shell:
    '''
    for f in !{stats_files} ; do
        sample="$(basename "${f}" .stats)"
        awk -v s="${sample}" 'BEGIN {OFS = "\t"} {print s, $0}' "${f}"
    done > !{basename}_per_sample_OTUs.stats
    '''
}
