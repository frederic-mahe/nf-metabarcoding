include { log_dir } from '../functions.nf'


process run_mumu {
    // [S43]: mumu (>=1.1.1) post-clustering filter. Inputs are the
    // reduced OTU table (amplicon + sample cols) and the self-search
    // match list; outputs are the new OTU table and the analysis log.
    //
    // [S45]: the mumu --log output is the canonical post-clustering
    // curation log. The intermediate _raw_mumu.table is **not**
    // published ([S46]); the publishDir pattern keeps the log only.
    // [D16]: logs go to logs/part_b/.
    publishDir path: { log_dir('part_b') }, mode: params.publish_mode, pattern: "*.log"

    input:
    path reduced_table
    path match_list
    val basename

    output:
    path "${reduced_table.baseName.replaceFirst(/_reduced$/, '_raw_mumu')}.table", emit: table
    path "${basename}_post_clustering_curation.log",                               emit: log

    shell:
    def new_table = "${reduced_table.baseName.replaceFirst(/_reduced$/, '_raw_mumu')}.table"
    // [S43]: mumu's relative-cooccurrence threshold is coupled to the
    // cleaving threshold --percentage ([S22]) as (1 - percentage), not an
    // independent knob: cleaving keeps a sub-seed present in >= percentage
    // of samples, and mumu merges a child OTU only when it co-occurs with
    // its parent in >= (1 - percentage) of the child's samples. The default
    // cleaving 0.05 gives 0.95. BigDecimal keeps the value exact (1 - 0.05
    // = 0.95) and free of binary-float noise.
    def minimum_relative_cooccurrence = BigDecimal.ONE - new BigDecimal(params.percentage.toString())
    """
    mumu \\
        --threads ${task.cpus} \\
        --otu_table ${reduced_table} \\
        --match_list ${match_list} \\
        --new_otu_table ${new_table} \\
        --minimum_relative_cooccurrence ${minimum_relative_cooccurrence} \\
        --log ${basename}_post_clustering_curation.log
    """

    stub:
    """
    touch ${reduced_table.baseName.replaceFirst(/_reduced$/, '_raw_mumu')}.table ${basename}_post_clustering_curation.log
    """
}
