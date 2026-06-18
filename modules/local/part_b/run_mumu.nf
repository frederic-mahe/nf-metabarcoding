include { normalize_path } from '../functions.nf'


process run_mumu {
    // [S43]: mumu (>=1.1.1) post-clustering filter. Inputs are the
    // reduced OTU table (amplicon + sample cols) and the self-search
    // match list; outputs are the new OTU table and the analysis log.
    //
    // [S45]: the mumu --log output is the canonical post-clustering
    // curation log. The intermediate _raw_mumu.table is **not**
    // published ([S46]); the publishDir pattern keeps the log only.
    publishDir path: { normalize_path(params.results_folder) }, mode: params.publish_mode, pattern: "*.log",
        enabled: params.results_folder != null

    input:
    path reduced_table
    path match_list
    val basename

    output:
    path "${reduced_table.baseName.replaceFirst(/_reduced$/, '_raw_mumu')}.table", emit: table
    path "${basename}_post_clustering_curation.log",                               emit: log

    shell:
    def new_table = "${reduced_table.baseName.replaceFirst(/_reduced$/, '_raw_mumu')}.table"
    """
    mumu \\
        --threads ${task.cpus} \\
        --otu_table ${reduced_table} \\
        --match_list ${match_list} \\
        --new_otu_table ${new_table} \\
        --log ${basename}_post_clustering_curation.log
    """
}
