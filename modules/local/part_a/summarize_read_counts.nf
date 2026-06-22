include { log_dir } from '../functions.nf'


process summarize_read_counts {
    // [S86]: port of the legacy genotoul per-sample read-tracking
    // summary. Gathers every regular sample's merging + trimming logs
    // and builds one TSV (samples / reads / assembled / F / R / passing
    // + a Total row). Published to logs/occurrence_table/ even though it
    // summarises Part A — it is a project-wide summary, one table per
    // run. bin/build_read_counts.sh reads the staged logs by filename.
    publishDir path: { log_dir('occurrence_table') }, mode: params.publish_mode, pattern: "*.tsv"

    input:
    path sample_ids        // sample_ids.txt, one regular sample id per line
    path logs              // every regular per-sample merging/trimming log
    val basename

    output:
    path "${basename}_read_counts.tsv", emit: table

    script:
    """
    build_read_counts.sh ${sample_ids} ${basename}_read_counts.tsv
    """

    stub:
    """
    touch ${basename}_read_counts.tsv
    """
}
