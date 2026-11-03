include { log_dir } from '../functions.nf'


process pool_fastq {
    // [S129]: pool the fastq files of a multi-file sample (folder
    // pooling, [S128], or repeated single-end samplesheet rows, [S132])
    // into one gzipped fastq file. part_A sorts the files by path; each
    // file is staged in its own numbered folder (in1/, in2/, ...) so
    // same-named files from different folders cannot collide, and the
    // list order is kept on the command line.
    publishDir path: { log_dir('part_a/per_sample') }, mode: params.publish_mode, pattern: "*.log"

    input:
    tuple val(sampleId), path(fastq_files, stageAs: 'in*/*')

    output:
    tuple val(sampleId), path("${sampleId}.fastq.gz"), emit: pooled
    path "${sampleId}_pooling.log",                    emit: log

    script:
    """
    pool_fastq.py \\
        --sample '${sampleId}' \\
        --output '${sampleId}.fastq.gz' \\
        --log '${sampleId}_pooling.log' \\
        ${fastq_files.collect { "'${it}'" }.join(' ')}
    """

    stub:
    """
    : | gzip > '${sampleId}.fastq.gz'
    touch '${sampleId}_pooling.log'
    """
}
