include { normalize_path } from '../functions.nf'


process compute_majority_assignment {
    // [S66]: opt-in final step of the regular Part C path. Recompute a
    // majority-rule taxonomy per OTU from the reference accessions in
    // the `references` column of the assigned table
    // (<basename>_table_assigned.tsv, [S51]) and the stampa-formatted
    // reference dataset ([S47]). Emits an independent three-column
    // table <basename>_table_assigned_majority.tsv
    // (OTU\tamplicon\ttaxonomy_majority). Never runs on the shadow
    // path; gated on params.majority_assignment and (by the startup
    // assert) on --taxonomy_method=stampa.
    publishDir path: { normalize_path(params.results_folder) }, mode: params.publish_mode,
        enabled: params.results_folder != null

    input:
    path assigned_table
    path reference_dataset
    val basename

    output:
    path "${basename}_table_assigned_majority.tsv"

    shell:
    '''
    majority_assignment.py \
        --input_table !{assigned_table} \
        --reference_db !{reference_dataset} \
        > !{basename}_table_assigned_majority.tsv
    '''
}
