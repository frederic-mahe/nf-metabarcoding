include { publish_dir } from '../functions.nf'


process extract_recluster_fasta {
    // [S105]/D-a: post-recluster sibling of `extract_mumu_fasta` ([S40])
    // — identical column layout (col 4 = amplicon, col 2 = abundance,
    // col 10 = sequence; header `>amplicon;size=total;`) and the same
    // `$2 != 0` safety filter. It runs only when the re-clustering pass
    // is on and re-extracts the published `<basename>_table.fas` from the
    // reclustered table, so the published table and FASTA stay
    // consistent. Kept a separate process (rather than an alias of
    // extract_mumu_fasta) so its publishDir is unconditional while
    // extract_mumu_fasta's is gated off on this path.
    publishDir path: { publish_dir('occurrence_table') }, mode: params.publish_mode

    input:
    path table

    output:
    path "${table.baseName}.fas", emit: fasta

    shell:
    '''
    awk 'NR > 1 && $2 != 0 {printf ">"$4";size="$2";\\n"$10"\\n"}' !{table} \
        > !{table.baseName}.fas
    '''

    stub:
    """
    touch ${table.baseName}.fas
    """
}
