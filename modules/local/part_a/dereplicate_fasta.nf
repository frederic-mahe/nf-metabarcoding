include { normalize_path } from '../functions.nf'


process dereplicate_fasta {
    // dereplicate and discard expected error values (ee)
    publishDir path: { normalize_path(params.fastq_folder) }, mode: params.publish_mode,
        enabled: params.fastq_folder != null

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
}
