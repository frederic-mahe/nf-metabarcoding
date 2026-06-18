include { normalize_path } from '../functions.nf'


process extract_mumu_fasta {
    // [S40]: post-mumu sibling of `extract_otu_fasta` — same column
    // layout, but skips rows whose `total == 0`. After the
    // size=0 → 1 awk hotfix in `rebuild_post_mumu_table` ([S44]) no
    // row carries `$2 == 0` anymore, so the filter is a no-op
    // safety net retained for byte parity with the legacy bash.
    publishDir path: { normalize_path(params.results_folder) }, mode: params.publish_mode,
        enabled: params.results_folder != null

    input:
    path table

    output:
    path "${table.baseName}.fas", emit: fasta

    shell:
    '''
    awk 'NR > 1 && $2 != 0 {printf ">"$4";size="$2";\\n"$10"\\n"}' !{table} \
        > !{table.baseName}.fas
    '''
}
