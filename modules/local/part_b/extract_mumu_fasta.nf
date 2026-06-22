include { publish_dir } from '../functions.nf'


process extract_mumu_fasta {
    // [S40]: post-mumu sibling of `extract_otu_fasta` — same column
    // layout, but skips rows whose `total == 0`. After the
    // size=0 → 1 awk hotfix in `rebuild_post_mumu_table` ([S44]) no
    // row carries `$2 == 0` anymore, so the filter is a no-op
    // safety net retained for byte parity with the legacy bash.
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
