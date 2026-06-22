include { publish_dir; log_dir } from '../functions.nf'


process dereplicate_fasta {
    // dereplicate and discard expected error values (ee)
    // [D15]: data (.fas) to per_sample/, log to logs/per_sample/.
    publishDir path: { publish_dir('per_sample') }, mode: params.publish_mode, pattern: "*.fas"
    publishDir path: { log_dir('per_sample') }, mode: params.publish_mode, pattern: "*.log"

    input:
    tuple val(sampleId), path(filtered_fasta)

    output:
    tuple val(sampleId), path("${sampleId}.fas"), emit: fasta
    path "${sampleId}_dereplicating.log",         emit: log

    shell:
    '''
    vsearch \
        --fastx_uniques !{filtered_fasta} \
        --sizeout \
        --quiet \
        --log !{sampleId}_dereplicating.log \
        --fasta_width 0 \
        --xee \
        --xlength \
        --fastaout - > !{sampleId}.fas
    '''

    stub:
    """
    touch ${sampleId}.fas ${sampleId}_dereplicating.log
    """
}
