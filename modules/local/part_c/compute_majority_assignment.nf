include { publish_dir } from '../functions.nf'


process compute_majority_assignment {
    // [S66]: opt-in final step of the regular Part C path. Recompute a
    // majority-rule taxonomy per OTU from the reference accessions in
    // the `references` column of the assigned table
    // (<basename>_table_assigned.tsv, [S51]) and the stampa-formatted
    // reference dataset ([S47]). Emits an independent three-column
    // table <basename>_taxonomy_stampa_majority.tsv
    // (OTU\tamplicon\ttaxonomy_majority). Never runs on the shadow
    // path; gated on params.majority_assignment and (by the startup
    // assert) on --taxonomy_method=stampa.
    publishDir path: { publish_dir('occurrence_table') }, mode: params.publish_mode

    input:
    path assigned_table
    path reference_dataset
    val basename

    output:
    path "${basename}_taxonomy_stampa_majority.tsv", emit: table

    shell:
    '''
    majority_assignment.py \
        --input_table !{assigned_table} \
        --reference_db !{reference_dataset} \
        > !{basename}_taxonomy_stampa_majority.tsv
    '''

    stub:
    """
    touch ${basename}_taxonomy_stampa_majority.tsv
    """
}
