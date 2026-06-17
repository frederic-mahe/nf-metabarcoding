include { normalize_path } from '../functions.nf'


process rebuild_post_mumu_table {
    // [S44]: wraps rebuild_table_after_mumu.py + the legacy
    // size=0 → 1 awk hotfix (so downstream `vsearch --sizein`
    // consumers don't choke on a zero-abundance row).
    //
    // [S46]: emits the final occurrence table as
    // <basename>_table.tsv.
    publishDir path: { normalize_path(params.results_folder) }, mode: params.publish_mode,
        enabled: params.results_folder != null

    input:
    path mumu_table
    path old_table
    val basename

    output:
    path "${basename}_table.tsv"

    shell:
    """
    rebuild_table_after_mumu.py \\
        --mumu_table ${mumu_table} \\
        --old_table  ${old_table} | \\
        awk 'BEGIN {FS = OFS = "\\t"} {if (\$2 == 0) {\$2 = 1} ; print \$0}' \\
        > ${basename}_table.tsv
    """
}
