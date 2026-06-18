include { normalize_path } from '../functions.nf'


process update_occurrence_table {
    // [S51]: splice the taxonomy assignments back onto Part B's
    // <basename>_table.tsv. The output is published as a sibling
    // file `<basename>_table_assigned.tsv` so Part B's unannotated
    // table is preserved alongside the annotated one (D04 sub-q2).
    publishDir path: { normalize_path(params.results_folder) }, mode: params.publish_mode,
        enabled: params.results_folder != null

    input:
    path occurrence_table
    path assignments
    val basename

    output:
    path "${basename}_table_assigned.tsv", emit: table

    shell:
    '''
    update_occurrence_table.py \
        --occurrence_table !{occurrence_table} \
        --assignments !{assignments} \
        > !{basename}_table_assigned.tsv
    '''
}
