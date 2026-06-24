include { publish_dir } from '../functions.nf'


process sort_taxonomy {
    // [S49]: stabilise the gathered stampa assignments. vsearch emits
    // each chunk's hits in a thread-dependent order and the per-chunk
    // slices are concatenated in a non-deterministic order, so the
    // merged file must be sorted to get reproducible output. The order
    // is the legacy stampa one: abundance descending, then amplicon
    // ascending (LC_ALL=C sort -k2,2nr -k1,1d on the
    // amplicon\tabundance\tidentity\ttaxonomy\treferences columns).
    //
    // This replaces the earlier collectFile(sort:) approach: that
    // operator's sort closure orders whole entries (chunks), not the
    // lines inside a chunk, so it only sorted correctly in the
    // degenerate one-record-per-chunk case and silently left
    // multi-record chunks (the default and the `local`/`demo` profiles)
    // unsorted. Delegating the sort to a real process fixes that and
    // also lifts the JVM-heap limit noted in the original Plan B.
    publishDir path: { publish_dir('occurrence_table') }, mode: params.publish_mode

    input:
    path taxonomy
    val basename

    output:
    path "${basename}_taxonomy_stampa.tsv", emit: taxonomy

    shell:
    '''
    LC_ALL=C sort -t "$(printf '\\t')" -k2,2nr -k1,1d !{taxonomy} \
        > !{basename}_taxonomy_stampa.tsv
    '''

    stub:
    """
    touch ${basename}_taxonomy_stampa.tsv
    """
}
